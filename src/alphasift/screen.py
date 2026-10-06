"""
Screen-level API: a directory of predictions in, a matrix and called sites out.

    from alphasift import analyze
    res = analyze("path/to/screen")
    res.M_pae        # partners x bait residues, raw min-PAE (A), NaN = no contact
    res.M_signal     # the same, transformed to (0, 1]; 0 = no contact
    res.sites        # called interaction sites, one dict per site

A "screen" is one bait predicted against many partners: the bait is whichever
chain sequence occurs in every prediction. The recurrence statistic needs that
design — it asks whether partners converge on the same residues more than
chance allows, which a single pair cannot answer. Single-pair input still
produces a matrix and strength regions, and `analyze` says so rather than
reporting a meaningless q-value.

Parameter defaults are the validated ones (see README); every one is exposed.
"""

from __future__ import annotations

import json
import tempfile
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from alphasift import readers
from alphasift.bait import (bait_from_chain, detect_bait,
                          detect_bait_from_partner_map, detect_bait_from_structures)
from alphasift.interfaces import merge_interfaces
from alphasift.readers import alphabridge as ab
from alphasift.readers.mmcif import chain_identities
from alphasift.matrix import bait_confidence_row, confidence_signal, contact_pairs
from alphasift.motifs import run_motif_engine
from alphasift.order import binder_clustering
from alphasift.partners import build_partner_annotations
from alphasift.rmsf import global_rmsf
from alphasift.stats import call_patches
from alphasift.structures import load_binder_structures

DEFAULTS = dict(
    # input
    source="pae",        # "pae" (any predictor) or "alphabridge" (AlphaFold3 only)
    bait_chain=None,     # name the bait chain instead of detecting it by sequence
    # matrix (alphasift.matrix)
    d0_contact=8.0,      # A; Cb-Cb bar that SELECTS candidate contacts
    d0_signal=10.0,      # A; fixed d0 of the pTM-style transform, all screens alike
    tau_strong_pae=10.0, # A; contact bar for the reported strength amount
    # statistics (alphasift.stats)
    P=5000, gap=1, min_size=4, pct=95, rng_seed=0,
    tau_single=0.8, min_single=4, single_gap=2,
    # where min_single applies: "partner" (>= 4 residues at >= tau_single anywhere on
    # the bait, then every sequence piece counts; 2026-09-30) or "group" (each piece)
    strength_scope="partner",
    # 3-D merge (alphasift.interfaces). "partner": sites LINK when they touch (support >=
    # merge_support), >= merge_min_shared partners bind both (signal >= partner_tau_in on
    # a site residue), shared/smaller >= merge_share_small and shared/larger >=
    # merge_share_large; interfaces = connected components. Short (< attach_min_len) or
    # weakly bound sites left alone group with each other and attach to one interface
    # (own partners at >= attach_tau_strong, else >= attach_tau_weak; partner-less sites
    # by closeness, support >= merge_prox_support). Any scipy linkage name (e.g.
    # "average", the rule before 2026-09-30) clusters on 1 - support instead.
    merge_3d=True, merge_method="partner", merge_support=0.50, d_merge=8.0,
    merge_min_shared=2, merge_share_small=0.7, merge_share_large=0.2, partner_tau_in=0.5,
    attach_tau_strong=0.8, attach_tau_weak=0.3, attach_min_len=4, merge_prox_support=0.9,
    # partner annotation (alphasift.partners): a partner is listed at an interface iff >= 1 pair
    # into it has s >= partner_tau_in (PAE <= 10 A); its residues are tiered H (s >=
    # partner_tau_high, PAE <= 5 A), M, L. Always from PAE, whichever `source` built the matrix.
    partner_tau_high=0.8,
    # motif engine (alphasift.motifs), on every interface: "c" (amino acids, BLOSUM62) or "3di_aa"
    # (amino acid + Foldseek 3Di of the partner chain at weight motif_w3di; needs Foldseek), or
    # None to skip. The validated viz_int settings (VALIDATION_LOG 2026-10-01/02).
    motif_engine="c", motif_n_perm=1000, motif_max_span=10, motif_touch_cb=7.0,
    motif_id_group=0.40, motif_w3di=1 / 3, motif_el_filter=True, foldseek_bin=None,
)


@dataclass
class Result:
    screen_dir: Path
    params: dict
    bait: dict
    preds: list                  # predictions contributing a row, in row order
    res_ids: np.ndarray          # bait residue numbers, 1..L
    M_pae: np.ndarray            # (n_preds, L) raw min-PAE, NaN = no contact
    M_signal: np.ndarray         # (n_preds, L) transformed, 0 = no contact
    source: str                  # which confidence source built the matrix
    sites: list                  # call_patches output, each with its 'interface'
    interfaces: object = None    # DataFrame: interface, segments, residues, n_patches,
                                 # channels, significant
    skipped: dict = field(default_factory=dict)   # prediction -> why no row
    single_pair: bool = False
    iface_res: dict = field(default_factory=dict)  # interface -> its bait residues
    partners: dict = field(default_factory=dict)   # interface -> [partner annotation] | None
    motifs: dict = field(default_factory=dict)     # motif engine output; {} when not run
    int_seq: dict = field(default_factory=dict)    # prediction -> {partner residue: 3-letter name}
    rmsf: tuple = (np.array([], int), np.array([], float))  # (bait residues, RMSF A) over binders
    identities: dict = field(default_factory=dict)  # prediction -> partner {uniprot, gene, ...} | None

    @property
    def s_sum(self):
        return self.M_signal.sum(0) if self.M_signal.size else np.zeros(len(self.res_ids))

    def row_order(self):
        """Leaf order for the matrix figure (weighted-Jaccard, average linkage)."""
        if self.M_signal.shape[0] < 3:
            return list(range(self.M_signal.shape[0])), None
        return binder_clustering(self.M_signal)

    def significant(self):
        return [s for s in self.sites if s["significant"]]


def _bait(screen: Path, bait_chain=None):
    """Bait detection, preferring what the screen actually provides.

    `bait_chain` names it outright and skips detection entirely. Otherwise: a
    `_partner_map.json` (database-sourced screens) names the bait; an
    AlphaBridge run carries chain sequences; failing both, sequences are read
    from the structures. The rule is the same in all three cases.
    """
    if bait_chain:
        return bait_from_chain(screen, bait_chain), "stated"
    partner_map = screen / "_partner_map.json"
    if partner_map.exists():
        return detect_bait_from_partner_map(screen, json.loads(partner_map.read_text())), "partner_map"
    if any(screen.glob("*/AlphaBridge/alphabridge_data.json")):
        try:
            return detect_bait(screen), "alphabridge"
        except (KeyError, IndexError, TypeError) as exc:
            # An AlphaBridge run is present but its JSON does not carry the chain
            # sequences this detector needs (older versions, or a partial run).
            # The structures are always there, so fall back rather than fail.
            warnings.warn(f"alphabridge_data.json has no usable chain sequences "
                          f"({exc!r}); identifying the bait from the structures instead")
    return detect_bait_from_structures(screen), "structures"


def load_contact_pairs(screen_dir, bait, preds, d0_contact=8.0):
    """{prediction: (bait_res, partner_res, pae_sym) | None} -- `contact_pairs` for each
    of `preds`, from the prediction files; where the raw PAE is not on disk, from the
    "contact_pairs" key of alphabridge_data.json (database-sourced screens). A prediction
    with neither maps to None, not to empty arrays: its partner side is UNKNOWN, not empty.
    Follows viz_int.confidence.load_screen_contact_pairs (2026-09-29)."""
    screen = Path(screen_dir)
    out = {}
    for pred in preds:
        p = readers.load_prediction(screen, pred)
        if p is None:
            abj = screen / pred / "AlphaBridge" / "alphabridge_data.json"
            cp = json.loads(abj.read_text()).get("contact_pairs") if abj.exists() else None
            if cp is None or cp.get("contact_d0") != d0_contact:
                out[pred] = None
            else:
                out[pred] = (np.asarray(cp["bait_res"], dtype=int),
                             np.asarray(cp["partner_res"], dtype=int),
                             np.asarray(cp["pae"], dtype=float))
            continue
        out[pred] = contact_pairs(p.chain_ids, p.res_ids, np.asarray(p.cb), np.asarray(p.pae),
                                  bait["chain_ids"][pred], d0_contact=d0_contact)
    return out


def load_matrix(screen_dir, d0_contact=8.0, d0_signal=10.0, bait_chain=None,
                source="pae"):
    """Return (bait, preds, M_pae, M_signal, res_ids, skipped).

    With `source="pae"` (the default, and the only option for AlphaFold2) the
    value is the symmetrised minimum PAE at each bait residue, taken per
    prediction from a precomputed `confidence_heatmap` in alphabridge_data.json
    where the raw PAE is no longer on disk, else from any registered reader.

    With `source="alphabridge"` it is instead the highest AlphaBridge sweep
    cut-off at which the residue belongs to an interface — available for
    AlphaFold3 screens that have been through an AlphaBridge run, since that
    sweep rests on AF3's distogram. `M_pae` is then all-NaN: that channel has no
    PAE of its own.

    Predictions with no readable value are skipped and reported in `skipped`,
    never silently dropped.
    """
    screen = Path(screen_dir)
    bait, _how = _bait(screen, bait_chain=bait_chain)
    chain_ids_map, L = bait["chain_ids"], bait["length"]

    if source not in ("pae", "alphabridge"):
        raise ValueError(f"source must be 'pae' or 'alphabridge', not {source!r}")

    preds, raw_rows, sig_rows, skipped = [], [], [], {}
    for pred, bc in chain_ids_map.items():
        if source == "alphabridge":
            j = ab.ab_json(screen, pred)
            if j is None:
                skipped[pred] = "no AlphaBridge output"
                continue
            preds.append(pred)
            raw_rows.append(np.full(L, np.nan))       # this channel has no PAE
            sig_rows.append(ab.bait_row(j, bc, L))
            continue
        raw = None
        ab_json = screen / pred / "AlphaBridge" / "alphabridge_data.json"
        if ab_json.exists():
            ch = json.loads(ab_json.read_text()).get("confidence_heatmap")
            if ch is not None:
                raw = np.array([np.nan if v is None else v for v in ch["raw_conf_pae"]],
                               dtype=float)
        if raw is None:
            p = readers.load_prediction(screen, pred)
            if p is None:
                skipped[pred] = "no readable structure/PAE pair"
                continue
            if p.pae.shape != (len(p.chain_ids), len(p.chain_ids)):
                skipped[pred] = (f"PAE {p.pae.shape} does not match {len(p.chain_ids)} tokens")
                continue
            bait_res_ids, raw_b, _interactor_chain = bait_confidence_row(
                p.chain_ids, p.res_ids, p.cb, p.pae, bc, d0_contact=d0_contact)
            raw = np.full(L, np.nan)
            for rid, v in zip(bait_res_ids, raw_b):
                if 1 <= rid <= L:
                    raw[rid - 1] = v
        preds.append(pred)
        raw_rows.append(raw)
        sig_rows.append(confidence_signal(raw, d0_signal))

    # Bait detection defines the screen from whichever files it read, so a
    # prediction it never saw would otherwise disappear without trace. Report
    # every prediction present on disk but absent from the matrix.
    for name in readers.discover(screen):
        if name not in preds and name not in skipped:
            skipped[name] = ("no AlphaBridge output" if source == "alphabridge"
                             else "not part of the detected bait screen")

    if source == "alphabridge" and skipped:
        # Recurrence tests convergence ACROSS partners, so a partial AlphaBridge
        # run does not merely lose rows -- it changes the test. Say so loudly.
        warnings.warn(
            f"source='alphabridge': {len(skipped)} of {len(skipped) + len(preds)} "
            "predictions have no AlphaBridge output and are absent from the matrix. "
            "Recurrence is a claim about convergence across partners, so it is being "
            "tested on a smaller screen than you have predictions for. Re-run "
            "AlphaBridge over the whole screen, or use source='pae'.")

    if source == "alphabridge" and not preds:
        raise FileNotFoundError(
            f"source='alphabridge' but no prediction under {screen} has "
            "AlphaBridge/alphabridge_data.json. That source needs an AlphaBridge run "
            "on AlphaFold3 output; use source='pae' otherwise.")

    M_pae = np.array(raw_rows) if raw_rows else np.zeros((0, L))
    M_signal = np.array(sig_rows) if sig_rows else np.zeros((0, L))
    return bait, preds, M_pae, M_signal, np.arange(1, L + 1), skipped


def analyze(screen_dir, **overrides) -> Result:
    """Build the matrix for one screen and call its interaction sites."""
    p = {**DEFAULTS, **overrides}
    screen_dir = Path(screen_dir)
    bait, preds, M_pae, M_signal, res_ids, skipped = load_matrix(
        screen_dir, d0_contact=p["d0_contact"], d0_signal=p["d0_signal"],
        bait_chain=p["bait_chain"], source=p["source"])

    if M_signal.shape[0] == 0:
        return Result(screen_dir, p, bait, preds, res_ids, M_pae, M_signal,
                      source=p["source"], sites=[], interfaces=_empty_interfaces(),
                      skipped=skipped)

    # Total confident engagement per bait residue, summed over partners: reported
    # as strength_A and never used for significance. It
    # is an Angstrom margin below the contact bar, so it exists only on the PAE
    # channel; the AlphaBridge sweep has no Angstrom scale to measure against.
    if p["source"] == "pae":
        ok = np.isfinite(M_pae) & (M_pae <= p["tau_strong_pae"])
        amount = np.where(ok, p["tau_strong_pae"] - np.nan_to_num(M_pae, nan=p["tau_strong_pae"]),
                          0.0).sum(0)
    else:
        amount = None

    sites = call_patches(M_signal, res_ids, P=p["P"], gap=p["gap"], min_size=p["min_size"],
                         pct=p["pct"], rng_seed=p["rng_seed"], strength_amount=amount,
                         tau_single=p["tau_single"],
                         min_single=p["min_single"], single_gap=p["single_gap"],
                         single_scope=p["strength_scope"])
    # ---- 3-D merge: which sites are one binding surface in space ------------
    binder_preds = sorted(pr for pr, ok in zip(preds, M_signal.sum(1) > 0) if ok)
    ifaces = _empty_interfaces()
    iface_res, partners, motifs, int_seq = {}, {}, {}, {}
    rmsf = (np.array([], int), np.array([], float))
    identities = {}
    if sites and binder_preds and p["merge_3d"] is not None:
        tbl = pd.DataFrame(sites)
        tbl["patch"] = tbl["res_start"].astype(str) + "-" + tbl["res_end"].astype(str)
        bait_atoms = {pr: readers.load_atoms(screen_dir, pr, bait["chain_ids"][pr])
                      for pr in binder_preds}
        # Which partners bind a site: each partner's matrix row, {bait residue: signal}
        # where it has a contact. viz_int rebuilds the same values from the contact
        # pairs (max of the transform over pairs = transform of the min PAE). With
        # source="alphabridge" the values are sweep cut-offs, so partner_tau_in and the
        # attach thresholds keep their numbers but not their PAE meaning.
        row = dict(zip(preds, M_signal))
        psig = {pr: {int(r): float(v) for r, v in zip(res_ids, row[pr]) if v > 0}
                for pr in binder_preds}
        ifaces, iface_res = merge_interfaces(
            tbl, bait_atoms, binder_preds, merge_3d=p["merge_3d"],
            merge_method=p["merge_method"], merge_support=p["merge_support"],
            d_merge=p["d_merge"], partner_sig=psig, min_shared=p["merge_min_shared"],
            share_small=p["merge_share_small"], share_large=p["merge_share_large"],
            tau_in=p["partner_tau_in"], tau_strong=p["attach_tau_strong"],
            tau_weak=p["attach_tau_weak"], min_len=p["attach_min_len"],
            prox_support=p["merge_prox_support"])
        where = {seg: row.interface for _, row in ifaces.iterrows() for seg in row.segments}
        for s in sites:
            s["interface"] = where.get((s["res_start"], s["res_end"]))

        # ---- partner side + motif engine (as viz_int.screen._run_structural_pipeline) ----
        # Both read the matrix's own contact pairs from PAE, whichever `source` built the matrix.
        pairs = load_contact_pairs(screen_dir, bait, binder_preds, d0_contact=p["d0_contact"])
        with tempfile.TemporaryDirectory() as tmp:
            cif_of = {pr: readers.model_cif(screen_dir, pr, tmp) for pr in binder_preds}
            want_cb = p["motif_engine"] is not None
            loaded = load_binder_structures(cif_of, bait, binder_preds, with_cb=want_cb)
            int_seq = loaded[2]
            partners = build_partner_annotations(ifaces, iface_res, pairs, int_seq,
                                                 d0_signal=p["d0_signal"], tau_in=p["partner_tau_in"],
                                                 tau_high=p["partner_tau_high"])
            # Every interface, significant or not (each result carries the flag): the bait-side
            # view exists whenever any partner touches the interface.
            if want_cb:
                eng_ifaces = [{"id": r["interface"], "residues": [int(x) for x in r["residues"]],
                               "significant": bool(r["significant"]),
                               "partners": [a["prediction"] for a in (partners.get(r["interface"]) or [])]}
                              for _, r in ifaces.iterrows()]
                if any(i["partners"] for i in eng_ifaces):
                    motifs = run_motif_engine(
                        p["motif_engine"], screen_dir.name, eng_ifaces, pairs, int_seq, cif_of, bait,
                        loaded[3], loaded[4], foldseek_bin=p["foldseek_bin"],
                        params=dict(N_PERM=p["motif_n_perm"], MAX_SPAN=p["motif_max_span"],
                                    TOUCH_CB=p["motif_touch_cb"], ID_GROUP=p["motif_id_group"],
                                    W3DI=p["motif_w3di"], EL_FILTER=p["motif_el_filter"]))
            rmsf = global_rmsf(cif_of, bait, binder_preds)
            identities = _identities(screen_dir, bait, cif_of)
    elif binder_preds:
        with tempfile.TemporaryDirectory() as tmp:
            cif_of = {pr: readers.model_cif(screen_dir, pr, tmp) for pr in binder_preds}
            rmsf = global_rmsf(cif_of, bait, binder_preds)
            identities = _identities(screen_dir, bait, cif_of)

    single = M_signal.shape[0] < 2
    if single:
        # Recurrence has no meaning with one partner: keep the regions, drop the claim.
        for s in sites:
            s["q_recur"], s["recurrent"], s["recurrent_core"] = None, False, []
            s["significant"] = s["strength_significant"]
    return Result(screen_dir, p, bait, preds, res_ids, M_pae, M_signal,
                  source=p["source"], sites=sites, interfaces=ifaces, skipped=skipped,
                  single_pair=single, iface_res=iface_res, partners=partners, motifs=motifs,
                  int_seq=int_seq, rmsf=rmsf, identities=identities)


def _identities(screen_dir, bait, cif_of):
    """Partner identity per prediction from its model header (mapped by sequence, see
    alphasift.readers.mmcif.chain_identities); the bait's own identity goes into
    `bait['identity']` from the first prediction that resolves it."""
    out = {}
    for pr, cif in cif_of.items():
        ids = chain_identities(Path(cif).read_text(), readers.load_sequences(screen_dir, pr) or {})
        bc = bait["chain_ids"][pr]
        out[pr] = next((v for ch, v in ids.items() if ch != bc), None)
        if bait.get("identity") is None and ids.get(bc):
            bait["identity"] = ids[bc]
    return out


def _empty_interfaces():
    return pd.DataFrame(columns=["interface", "n_patches", "segments", "residues",
                                 "channels", "significant"])
