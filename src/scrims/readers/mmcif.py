"""
mmCIF coordinates (Cbeta per residue) and AlphaFold3 confidences.json (PAE).

Extracted verbatim from viz_int.contact_proxy (`pseudo_cb`, `Residues`,
`parse_cif_cb`, `Confidences`, `load_confidences`). The contact-proxy machinery
that shared that module -- the legacy geometry channel -- is out of scope.
Self-contained: no Bio.PDB, no AlphaBridge.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def pseudo_cb(n: np.ndarray, ca: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Ideal Cbeta position from backbone N, CA, C (used for glycine).

    Standard tetrahedral placement. If you have a preferred constant set,
    replace this function body; nothing else depends on the constants.
    """
    b = ca - n
    v = c - ca
    a = np.cross(b, v)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * v + ca


# --------------------------------------------------------------------------- #
# mmCIF coordinate parsing (Cbeta per residue)
# --------------------------------------------------------------------------- #

@dataclass
class Residues:
    """Per-residue Cbeta coordinates, keyed for alignment with PAE tokens."""

    chain_ids: np.ndarray   # (N,) label_asym_id, str
    res_ids: np.ndarray     # (N,) label_seq_id, int
    res_names: np.ndarray   # (N,) label_comp_id, str
    cb: np.ndarray          # (N, 3) float, Cbeta (pseudo for Gly)

    def index(self) -> dict[tuple[str, int], int]:
        return {
            (str(cid), int(rid)): i
            for i, (cid, rid) in enumerate(zip(self.chain_ids, self.res_ids))
        }


def _atom_site_columns(lines: list[str], start: int) -> tuple[dict[str, int], int]:
    """Read a loop_ header starting at `start`, return {column_name: index}
    and the line index of the first data row."""
    cols: dict[str, int] = {}
    i = start
    while i < len(lines) and lines[i].strip().startswith("_atom_site."):
        name = lines[i].strip().split(".", 1)[1]
        cols[name] = len(cols)
        i += 1
    return cols, i


def parse_cif_cb(cif_path: str | Path) -> Residues:
    """Parse an mmCIF file and return one Cbeta coordinate per polymer residue.

    Cbeta is the real CB atom when present, otherwise a pseudo-Cbeta computed
    from N, CA, C (glycine). Residues lacking the backbone atoms needed for a
    pseudo-Cbeta are dropped with their coordinate set to NaN.
    """
    lines = Path(cif_path).read_text().splitlines()

    # Locate the _atom_site loop header.
    cols = None
    data_start = None
    for i, ln in enumerate(lines):
        if ln.strip() == "loop_":
            j = i + 1
            # skip blank lines
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines) and lines[j].strip().startswith("_atom_site."):
                cols, data_start = _atom_site_columns(lines, j)
                break
    if cols is None:
        raise ValueError(f"No _atom_site loop found in {cif_path}")

    required = [
        "label_asym_id", "label_seq_id", "label_comp_id",
        "label_atom_id", "Cartn_x", "Cartn_y", "Cartn_z",
    ]
    missing = [c for c in required if c not in cols]
    if missing:
        raise ValueError(f"_atom_site loop missing columns {missing} in {cif_path}")

    ic = cols["label_asym_id"]
    ir = cols["label_seq_id"]
    icomp = cols["label_comp_id"]
    iatom = cols["label_atom_id"]
    ix, iy, iz = cols["Cartn_x"], cols["Cartn_y"], cols["Cartn_z"]
    ncol = len(cols)

    # Gather backbone/CB atoms per (chain, res_id), preserving first-seen order.
    order: list[tuple[str, int]] = []
    atoms: dict[tuple[str, int], dict] = {}

    for ln in lines[data_start:]:
        s = ln.strip()
        if not s or s.startswith("#") or s.startswith("_") or s == "loop_":
            if s.startswith("_atom_site."):
                continue
            # a new category / loop begins -> atom_site block is done
            if s and not s[0].isspace() and (s.startswith("_") or s == "loop_" or s == "#"):
                break
            continue
        f = s.split()
        if len(f) < ncol:
            continue
        if f[0] not in ("ATOM", "HETATM"):
            continue
        chain = f[ic]
        raw_res = f[ir]
        if raw_res in (".", "?"):
            continue  # non-polymer token without a sequence id
        res_id = int(raw_res)
        atom = f[iatom].strip('"')
        key = (chain, res_id)
        if key not in atoms:
            atoms[key] = {"comp": f[icomp]}
            order.append(key)
        if atom in ("N", "CA", "C", "CB"):
            atoms[key][atom] = np.array(
                [float(f[ix]), float(f[iy]), float(f[iz])], dtype=float
            )

    chain_ids, res_ids, res_names, cbs = [], [], [], []
    for key in order:
        chain, res_id = key
        rec = atoms[key]
        if "CB" in rec:
            cb = rec["CB"]
        elif all(k in rec for k in ("N", "CA", "C")):
            cb = pseudo_cb(rec["N"], rec["CA"], rec["C"])
        else:
            cb = np.array([np.nan, np.nan, np.nan])
        chain_ids.append(chain)
        res_ids.append(res_id)
        res_names.append(rec["comp"])
        cbs.append(cb)

    return Residues(
        chain_ids=np.array(chain_ids),
        res_ids=np.array(res_ids, dtype=int),
        res_names=np.array(res_names),
        cb=np.vstack(cbs) if cbs else np.zeros((0, 3)),
    )


# --------------------------------------------------------------------------- #
# PAE / native-contact loading
# --------------------------------------------------------------------------- #

@dataclass
class Confidences:
    pae: np.ndarray                 # (T, T)
    token_chain_ids: np.ndarray     # (T,)
    token_res_ids: np.ndarray       # (T,)
    contact_probs: np.ndarray | None  # (T, T) native ground truth, if present


def load_confidences(
    confidences_json: str | Path,
    *,
    pae_key: str = "pae",
    contact_key: str = "contact_probs",
    chain_key: str = "token_chain_ids",
    res_key: str = "token_res_ids",
) -> Confidences:
    """Load PAE (and native contact_probs if present) from an AF3
    confidences.json. Key names are overridable in case a variant differs."""
    d = json.loads(Path(confidences_json).read_text())
    for k in (pae_key, chain_key, res_key):
        if k not in d:
            raise KeyError(
                f"'{k}' not found in {confidences_json}. "
                f"Available keys: {sorted(d)}"
            )
    pae = np.asarray(d[pae_key], dtype=float)
    contact = np.asarray(d[contact_key], dtype=float) if contact_key in d else None
    return Confidences(
        pae=pae,
        token_chain_ids=np.asarray(d[chain_key]),
        token_res_ids=np.asarray(d[res_key], dtype=int),
        contact_probs=contact,
    )



def parse_cif_atoms(cif_path: str | Path, chain: str | None = None,
                    text: str | None = None) -> dict:
    """{residue number: (n_atoms, 3) heavy-atom coordinates} for one chain.

    The 3-D interface merge needs heavy atoms, not just Cbeta. Written against
    the same `_atom_site` loop `parse_cif_cb` reads, so SCRIMS needs no second
    structure library; `chain=None` returns {chain: {res: coords}} for all chains.

    Only `ATOM` records are kept, which is how the reference implementation's
    Bio.PDB pass behaved (it filtered on a blank het-flag). Modified residues
    written as HETATM are therefore dropped here too — none occur in the
    AlphaFold output this was validated on, but a screen containing PTMs would
    lose those residues from the merge and should be checked.
    """
    lines = (text if text is not None else Path(cif_path).read_text()).splitlines()
    cols, start = None, None
    for i, ln in enumerate(lines):
        if ln.strip().startswith("_atom_site."):
            cols, j = {}, i
            while j < len(lines) and lines[j].strip().startswith("_atom_site."):
                cols[lines[j].strip().split(".", 1)[1]] = len(cols)
                j += 1
            start = j
            break
    if cols is None:
        raise ValueError(f"no _atom_site loop in {cif_path}")

    ig = cols.get("group_PDB")
    ic, ir = cols["label_asym_id"], cols["label_seq_id"]
    ia = cols["label_atom_id"]
    ix, iy, iz = cols["Cartn_x"], cols["Cartn_y"], cols["Cartn_z"]
    ncol = len(cols)

    out: dict[str, dict[int, list]] = {}
    for ln in lines[start:]:
        s = ln.strip()
        if not s or s[0] == "#" or s.startswith("_") or s == "loop_":
            if s.startswith(("data_", "loop_", "_")) or not s:
                continue
            break
        f = s.split()
        if len(f) < ncol:
            continue
        if ig is not None and f[ig] != "ATOM":
            continue
        ch = f[ic]
        if chain is not None and ch != chain:
            continue
        rid = f[ir]
        if rid in (".", "?"):
            continue
        atom = f[ia].strip('"')
        if atom.startswith("H"):
            continue
        out.setdefault(ch, {}).setdefault(int(rid), []).append(
            (float(f[ix]), float(f[iy]), float(f[iz])))

    packed = {c: {r: np.asarray(v, dtype=float) for r, v in d.items()} for c, d in out.items()}
    if chain is not None:
        return packed.get(chain, {})
    return packed


def _header_loops(text: str, prefixes=("_entity_poly.", "_struct_ref.")) -> dict:
    """{category: [row dicts]} for the header loops named by `prefixes` (quoted values
    with spaces handled). Only multi-row `loop_` blocks are read."""
    import shlex
    out, lines, i = {}, text.splitlines(), 0
    while i < len(lines):
        if lines[i].strip() == "loop_" and i + 1 < len(lines) and lines[i + 1].startswith(prefixes):
            cat = lines[i + 1].split(".")[0]
            cols = []
            i += 1
            while i < len(lines) and lines[i].startswith(cat + "."):
                cols.append(lines[i].split(".", 1)[1].strip())
                i += 1
            rows = []
            while i < len(lines) and lines[i].strip() not in ("#", "loop_", "") and not lines[i].startswith("_"):
                rows.append(shlex.split(lines[i]))
                i += 1
            out[cat] = [dict(zip(cols, r)) for r in rows]
            continue
        i += 1
    return out


def chain_identities(text: str, chain_sequences: dict) -> dict:
    """{chain: {'uniprot', 'entry_name', 'gene', 'organism', 'description'} | None} from an
    mmCIF header (`gene`/`organism` from ModelCIF's `_ma_target_ref_db_details`, when present).

    Maps each chain BY SEQUENCE: its full sequence must equal exactly one
    `_entity_poly` sequence, whose entity is then looked up in `_struct_ref` and
    `_entity`. The header's chain -> entity link (`_struct_asym`) is never used: in the
    raw NVIDIA files it is wrong in ~37 % of models (viz_int docs/PROTEOME_SCALING.md).
    Checked: the bait chain resolves to the bait's own accession in 300/300 models
    (RHOQ, VAPA). Chains that do not resolve map to None; AlphaFold3 server output
    carries no `_struct_ref`, so everything is None there.
    """
    h = _header_loops(text, ("_entity_poly.", "_struct_ref.", "_entity.", "_ma_target_ref_db_details."))
    seq_of = {r.get("entity_id"): r.get("pdbx_seq_one_letter_code", "").replace("\n", "")
              for r in h.get("_entity_poly", [])}
    ref = {r.get("entity_id"): r for r in h.get("_struct_ref", [])}
    desc = {r.get("id"): r.get("pdbx_description") for r in h.get("_entity", [])}
    ma = {r.get("target_entity_id"): r for r in h.get("_ma_target_ref_db_details", [])}
    out = {}
    for ch, s in chain_sequences.items():
        ents = [e for e, v in seq_of.items() if v == s]
        r = ref.get(ents[0]) if len(ents) == 1 else None
        m = ma.get(ents[0], {}) if r is not None else {}
        out[ch] = None if r is None else {"uniprot": r.get("pdbx_db_accession"),
                                          "entry_name": r.get("db_code"),
                                          "gene": m.get("gene_name"),
                                          "organism": m.get("organism_scientific"),
                                          "description": desc.get(ents[0])}
    return out
