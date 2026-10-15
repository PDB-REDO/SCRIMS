"""
Machine-readable output: one JSON document per screen (`scrims.json`) holding
everything the figures and the HTML report are drawn from, plus flat TSV tables.

Schema "scrims/1" (top-level keys):
  schema, scrims_version
  screen        {name, dir, source, single_pair, n_predictions, skipped: {prediction: reason}}
  bait          {sequence, length, identity, chain_ids: {prediction: chain}}
  params        every parameter the run used
  predictions   [{id, row, binder, identity}]  -- `row` indexes the matrix; identity =
                {uniprot, entry_name, gene, organism, description} from the model header, or null
  matrix        {residues, row_order, pae, signal}: pae = raw min-PAE in A (2 decimals,
                null = no contact; all null on the AlphaBridge source), signal =
                1/(1+(PAE/d0)^2) in (0, 1], 0 = no contact (4 decimals); row_order = the
                display order of the BINDER rows (>= 1 contact), clustered among themselves
                as on the website (indices into rows; non-binders are not listed)
  tracks        {summed_signal: per residue (what recurrence tests), rmsf: {residues, values}}
  sites         one per called site, as sites.tsv
  interfaces    [{id, significant, segments, residues, n_patches, channels,
                  partners: [{id, interactor, identity, ranges: [[start, end, tier, letters]],
                              n_high, n_medium, n_low, bait_residues, n_pairs, best_pae,
                              median_pae}] | null (no pairwise PAE: unknown, not empty),
                  motif: null (no partner touches it, or the engine did not run) |
                         {panel1: [{pos, n_groups, agree, q, letters, letters_grouped}],
                          (letters = share of partner residues at the bait residue; letters_grouped =
                           the same with each partner group of >= motif_id_group identity weighing 1
                           in total, split over its members -- the bait-side classic logo's input)
                          elements: [{positions, n_groups, q, letters, orientation, consistency,
                                      carriers: [{id, group, residues, map}]}],
                          motifs: [{id, elements, positions, n_groups, carriers, bipartite,
                                    linked_by, touching, candidates: [{elm, name, tier,
                                    groups_matched, groups, expected, p, domain, matches}]}],
                          (candidates = known ELM classes its partners match beyond chance,
                           scrims.motif_candidates; [] = none suggested)
                          motif_relations: [{a, b, relation, shared_positions, k, exclusive,
                                             q_dep, status}],
                          bipartite, alternative, control: {false_elements, false_columns},
                          partners_used, partners_total,
                          alignments: {motif id: scrims.motif_align.align_motif(...)}}}]
  motif_engine  {engine, params, groups: {prediction: group id}} | null

The per-interface `motif` block has the shape of the companion website's
`motif_elements_confidence` entries (viz_int webexport, schema 0.9), without its
website-only fields, so the same reading applies to both.
"""

from __future__ import annotations

import csv
import json
import math

import numpy as np

from scrims.motif_candidates import motif_candidates, parse_domains

SCHEMA = "scrims/1"


def _f(x, nd=None):
    if x is None:
        return None
    x = float(x)
    if not math.isfinite(x):
        return None
    return round(x, nd) if nd is not None else x


def _touching(touch, m, S2):
    """{carrier: {bait position: sorted partner residues touching it}} over motif m's bait positions, for every carrier
    of m's elements (viz_int 2026-10-08: keyed by bait position for the binding sub-groups)."""
    pos = {int(x) for x in m["positions"]}
    cars = sorted({pr for i in m["elements"] for pr in S2["elements"][i]["residues"]} | set(m["carriers"]))
    return {c: {str(b): sorted(int(r) for r in rs) for b, rs in sorted((touch.get(c) or {}).items()) if int(b) in pos and rs}
            for c in cars}


def _motif_block(r, groups, seed=None, int_seq=None, bait_domains=frozenset()):
    """One interface's engine result in the website's 0.9 entry shape (viz_int
    webexport._motif_elements_confidence, minus gene/uniprot/cif per carrier)."""
    S2 = r["panel2"]["shuffle"]
    els = []
    for el in S2["elements"]:
        els.append({"positions": [int(x) for x in el["positions"]], "n_groups": int(el["n_groups"]),
                    "q": _f(el["q"]), "letters": el["letters"], "orientation": el["orientation"],
                    "consistency": _f(el["consistency"]),
                    "carriers": [{"id": pr, "group": groups.get(pr, pr),
                                  "residues": [int(x) for x in el["residues"][pr]],
                                  "map": {str(k): int(v) for k, v in el["map"][pr].items()}}
                                 for pr in sorted(el["residues"])]})

    def pair_rows(rows, qk):
        return [{"a": int(x["a"]), "b": int(x["b"]), "k": int(x["k"]), "q": _f(x[qk]),
                 "linkers": [[o, int(n)] for o, n in x["linkers"]],
                 **({"fixed_linker": bool(x["fixed_linker"])} if "fixed_linker" in x else {})}
                for x in rows]
    return {
        "panel1": [{"pos": int(c["pos"]), "n_groups": int(c["n_groups"]), "agree": _f(c["agree"]),
                    "q": _f((c["q"] or {}).get("shuffle")), "letters": c["letters"],
                    "letters_grouped": c.get("letters_grouped", {})} for c in r["panel1"]],
        "elements": els,
        "motifs": [{"id": m["id"], "elements": [int(x) for x in m["elements"]],
                    "positions": [int(x) for x in m["positions"]], "n_groups": int(m["n_groups"]),
                    "carriers": list(m["carriers"]), "bipartite": bool(m["bipartite"]),
                    "linked_by": list(m.get("linked_by", [])),
                    # every carrier residue touching the motif's bait positions (display only; viz_int 2026-10-07)
                    **({"touching": _touching(r["touch"], m, S2)} if "touch" in r else {}),
                    # known ELM class its partners match beyond chance (viz_int IMPROVEMENTS #41; 2026-10-12)
                    **({"candidates": motif_candidates(seed, r["_iface"], m, S2["elements"], r.get("touch") or {},
                                                       int_seq, groups, bait_domains)} if int_seq else {})}
                   for m in S2.get("motifs", [])],
        "motif_relations": [{"a": x["a"], "b": x["b"], "relation": x["relation"],
                             "shared_positions": [int(y) for y in x["shared_positions"]],
                             "k": int(x["k"]), "exclusive": bool(x["exclusive"]), "q_dep": _f(x["q_dep"]),
                             "status": x.get("status", "unresolved")} for x in S2.get("motif_relations", [])],
        "bipartite": pair_rows(S2["bipartite"], "q_enr"),
        "alternative": pair_rows(S2["alternative"], "q_dep"),
        "control": {"false_elements": int(r["control"]["elements"]),
                    "false_columns": [int(x) for x in r["control_panel1"]]},
        "partners_used": int(r["partners_used"]), "partners_total": int(r["partners_total"]),
    }


def document(res):
    """The `scrims/1` document for one analysed screen (an scrims.Result)."""
    from scrims import __version__
    from scrims.motif_align import align_motif

    # Display order as on the website: binder rows only (>= 1 contact), clustered among themselves.
    from scrims.order import binder_clustering
    b_idx = [k for k, ok in enumerate(res.M_signal.sum(1) > 0) if ok] if res.M_signal.size else []
    row_order = [b_idx[k] for k in binder_clustering(res.M_signal[b_idx])[0]] if b_idx else []
    eng = res.motifs or {}
    groups = eng.get("groups", {})
    seed = (res.bait.get("identity") or {}).get("gene") or res.screen_dir.name
    bait_domains = parse_domains(res.params.get("bait_domains", ""))
    binders = set(res.identities) | {pr for pr, ok in zip(res.preds, res.M_signal.sum(1) > 0) if ok} \
        if res.M_signal.size else set()

    ifaces = []
    for _, row in res.interfaces.iterrows():
        iid = row["interface"]
        anns = res.partners.get(iid) if res.partners else None
        partners = None if anns is None else [
            {"id": a["prediction"], "interactor": a["interactor"], "identity": res.identities.get(a["prediction"]),
             "ranges": a["ranges"], "n_high": a["n_high"], "n_medium": a["n_medium"], "n_low": a["n_low"],
             "bait_residues": a["bait_residues"], "n_pairs": a["n_pairs"],
             "best_pae": _f(a["best_pae"]), "median_pae": _f(a["median_pae"])} for a in anns]
        r = (eng.get("ifaces") or {}).get(iid)
        motif = None
        if r is not None:
            motif = _motif_block({**r, "_iface": iid}, groups, seed, res.int_seq, bait_domains)
            motif["alignments"] = {m["id"]: align_motif(motif, m, partners) for m in motif["motifs"]}
        ifaces.append({"id": iid, "significant": bool(row["significant"]),
                       "segments": [[int(a), int(b)] for a, b in row["segments"]],
                       "residues": [int(x) for x in row["residues"]], "n_patches": int(row["n_patches"]),
                       "channels": row["channels"], "partners": partners, "motif": motif})

    pae = [[None if not np.isfinite(v) else round(float(v), 2) for v in r] for r in res.M_pae]
    sig = [[round(float(v), 4) for v in r] for r in res.M_signal]
    g_res, g_rmsf = res.rmsf
    return {
        "schema": SCHEMA, "scrims_version": __version__,
        "screen": {"name": res.screen_dir.name, "dir": str(res.screen_dir), "source": res.source,
                   "single_pair": res.single_pair, "n_predictions": len(res.preds), "skipped": res.skipped},
        "bait": {"sequence": res.bait.get("sequence"), "length": int(res.bait["length"]),
                 "identity": res.bait.get("identity"), "chain_ids": res.bait.get("chain_ids", {})},
        "params": {k: (v if not isinstance(v, float) or math.isfinite(v) else None) for k, v in res.params.items()},
        "predictions": [{"id": pr, "row": i, "binder": pr in binders, "identity": res.identities.get(pr)}
                        for i, pr in enumerate(res.preds)],
        "matrix": {"residues": [int(x) for x in res.res_ids], "row_order": [int(x) for x in row_order],
                   "pae": pae, "signal": sig},
        "tracks": {"summed_signal": [round(float(x), 4) for x in res.s_sum],
                   "rmsf": {"residues": [int(x) for x in g_res], "values": [round(float(x), 3) for x in g_rmsf]}},
        "sites": [{k: (_f(v) if isinstance(v, float) else v) for k, v in s.items()} for s in res.sites],
        "interfaces": ifaces,
        "motif_engine": ({"engine": eng.get("engine"), "params": eng.get("params"), "groups": groups}
                         if eng else None),
    }


def label(identity, pred):
    """Display name for a partner: its gene, else the UniProt entry name, else the prediction id."""
    if identity:
        return identity.get("gene") or (identity.get("entry_name") or "").split("_")[0] or pred
    return pred


def write_tables(doc, outdir):
    """partners.tsv, bait_side.tsv and motifs.tsv from an `scrims/1` document."""
    def tsv(name, header, rows):
        with open(outdir / name, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(header)
            w.writerows(rows)
    rows = []
    for i in doc["interfaces"]:
        for a in i["partners"] or []:
            idn = a["identity"] or {}
            rows.append([i["id"], a["id"], label(a["identity"], a["id"]), idn.get("uniprot", ""),
                         a["n_high"], a["n_medium"], a["n_low"], a["best_pae"], a["median_pae"],
                         ";".join(f"{s}-{e}:{t}:{q}" for s, e, t, q in a["ranges"]),
                         ";".join(map(str, a["bait_residues"]))])
    tsv("partners.tsv", ["interface", "prediction", "partner", "uniprot", "n_high", "n_medium", "n_low",
                         "best_pae", "median_pae", "ranges(start-end:tier:sequence)", "bait_residues"], rows)
    rows = []
    for i in doc["interfaces"]:
        for c in (i["motif"] or {}).get("panel1", []):
            rows.append([i["id"], c["pos"], doc["bait"]["sequence"][c["pos"] - 1] if doc["bait"]["sequence"] else "",
                         c["n_groups"], "" if c["agree"] is None else c["agree"],
                         "" if c["q"] is None else c["q"],
                         ";".join(f"{a}:{v}" for a, v in sorted(c["letters"].items(), key=lambda kv: -kv[1])),
                         ";".join(f"{a}:{v}" for a, v in sorted(c.get("letters_grouped", {}).items(), key=lambda kv: -kv[1]))])
    tsv("bait_side.tsv", ["interface", "bait_residue", "bait_aa", "n_groups", "agreement", "q", "partner_letters",
                          "partner_letters_grouped"], rows)
    rows = []
    for i in doc["interfaces"]:
        m = i["motif"] or {}
        for mo in m.get("motifs", []):
            al = m["alignments"][mo["id"]]
            names = [label(next((p["identity"] for p in i["partners"] or [] if p["id"] == c), None), c)
                     for c in mo["carriers"]]
            # 2026-10-09: the motif consensus (motif_align.motif_consensus; before: the top letter of every column)
            cons = "-".join(t for _, t in al.get("consensus") or [])
            rows.append([i["id"], i["significant"], mo["id"], ";".join(map(str, mo["positions"])), mo["n_groups"],
                         len(mo["carriers"]), ",".join(mo["linked_by"]), mo["bipartite"], cons,
                         al.get("quality", ""), al.get("quality_top", ""),
                         min((m["elements"][k]["q"] for k in mo["elements"] if m["elements"][k]["q"] is not None), default=""),
                         "antiparallel" if al.get("flipped") else "",
                         ";".join(f"{c['elm']} ({c['tier']}, {c['groups_matched']}/{c['groups']} groups, p={c['p']:.2g}, "
                                  f"domain {c['domain']})" for c in mo.get("candidates") or []),
                         ";".join(names)])
    tsv("motifs.tsv", ["interface", "interface_significant", "motif", "bait_residues", "n_groups", "n_carriers",
                       "linked_by", "bipartite", "partner_consensus", "quality_bits", "quality_top_pct", "motif_q", "orientation",
                       "known_class", "carriers"], rows)


def write_json(doc, path):
    path.write_text(json.dumps(doc, separators=(",", ":")))
