"""
AlphaFold2 output in AFDB/ColabFold bulk packaging: `<id>-model_v1.cif[.zst]`
plus `<id>-predicted_aligned_error_v1.json[.zst]`.

Extracted from viz_int.colabfold_adapter (`unzst`, `parse_cif`, `parse_pae`,
`load_entry`). The single-model geometry contact matrix that module also
provided is out of scope. zstandard is only imported when a `.zst` file is
actually opened, so plain (uncompressed) AFDB downloads work without it.
"""

from __future__ import annotations

import json
from dataclasses import dataclass

import numpy as np


def read_maybe_zst(path: str) -> bytes:
    """File contents, transparently decompressing `.zst`. Plain AFDB downloads
    therefore need no zstandard install."""
    if str(path).endswith(".zst"):
        import zstandard
        with open(path, "rb") as f:
            return zstandard.ZstdDecompressor().stream_reader(f).read()
    with open(path, "rb") as f:
        return f.read()


def _pseudo_cb(n, ca, c):
    b = ca - n
    v = c - ca
    a = np.cross(b, v)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * v + ca


@dataclass
class Entry:
    token_chain_ids: np.ndarray   # (T,)
    token_res_ids: np.ndarray     # (T,)
    comp_ids: np.ndarray          # (T,) 3-letter residue names
    cb: np.ndarray                # (T, 3) Cβ (pseudo for Gly)
    plddt: np.ndarray             # (T,) per-residue mean pLDDT (0-100)
    pae: np.ndarray               # (T, T)


def parse_cif(cif_text: str):
    """Parse the mmCIF _atom_site loop → per-residue Cβ, pLDDT, chain/res ids
    in token order (as they appear, i.e. chain by chain)."""
    lines = cif_text.splitlines()
    # locate the _atom_site loop header
    cols, start = None, None
    for i, ln in enumerate(lines):
        if ln.strip() == "loop_":
            j = i + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines) and lines[j].strip().startswith("_atom_site."):
                cols = {}
                while j < len(lines) and lines[j].strip().startswith("_atom_site."):
                    cols[lines[j].strip().split(".", 1)[1]] = len(cols)
                    j += 1
                start = j
                break
    if cols is None:
        raise ValueError("no _atom_site loop")

    ic, ir, icomp, iatom = (cols["label_asym_id"], cols["label_seq_id"],
                            cols["label_comp_id"], cols["label_atom_id"])
    ix, iy, iz = cols["Cartn_x"], cols["Cartn_y"], cols["Cartn_z"]
    ib = cols["B_iso_or_equiv"]
    ncol = len(cols)

    order, atoms = [], {}
    for ln in lines[start:]:
        s = ln.strip()
        if not s or s[0] == "#" or s.startswith("_") or s == "loop_":
            if s and not s.startswith("_atom_site."):
                break
            continue
        f = s.split()
        if len(f) < ncol or f[0] not in ("ATOM", "HETATM"):
            continue
        rid = f[ir]
        if rid in (".", "?"):
            continue
        key = (f[ic], int(rid))
        atom = f[iatom].strip('"')
        if key not in atoms:
            atoms[key] = {"comp": f[icomp], "b": [], "xyz": {}}
            order.append(key)
        atoms[key]["b"].append(float(f[ib]))
        if atom in ("N", "CA", "C", "CB"):
            atoms[key]["xyz"][atom] = np.array([float(f[ix]), float(f[iy]), float(f[iz])])

    chain_ids, res_ids, comps, cbs, plddts = [], [], [], [], []
    for key in order:
        rec = atoms[key]
        xyz = rec["xyz"]
        if "CB" in xyz:
            cb = xyz["CB"]
        elif all(k in xyz for k in ("N", "CA", "C")):
            cb = _pseudo_cb(xyz["N"], xyz["CA"], xyz["C"])
        else:
            cb = np.array([np.nan, np.nan, np.nan])
        chain_ids.append(key[0]); res_ids.append(key[1]); comps.append(rec["comp"])
        cbs.append(cb); plddts.append(float(np.mean(rec["b"])))
    return (np.array(chain_ids), np.array(res_ids, dtype=int), np.array(comps),
            np.vstack(cbs), np.array(plddts))


def parse_pae(pae_json: bytes) -> np.ndarray:
    d = json.loads(pae_json)
    if isinstance(d, list):
        d = d[0]
    return np.asarray(d["predicted_aligned_error"], dtype=float)


def load_entry(cif_zst_path: str, pae_zst_path: str) -> Entry:
    chain_ids, res_ids, comps, cb, plddt = parse_cif(read_maybe_zst(cif_zst_path).decode())
    pae = parse_pae(read_maybe_zst(pae_zst_path))
    T = len(chain_ids)
    if pae.shape != (T, T):
        raise ValueError(f"PAE {pae.shape} != ({T},{T}); CIF residues and PAE tokens misaligned")
    return Entry(chain_ids, res_ids, comps, cb, plddt, pae)

