"""
Per-residue flexibility of the bait across the screen: all-heavy-atom RMSF over
the binders' bait chains after a reference-free consensus superposition
(Theseus-style). Drawn as the RMSF track above the matrix, as on the website.

`_kabsch_Rt`, `_consensus` and `global_rmsf` are extracted verbatim from
viz_int.plotdata (ov_4), except that each prediction's model file is passed in
(`cif_of`, see alphasift.readers.model_cif) instead of globbed from the screen.
"""

from __future__ import annotations

import numpy as np
from numpy.linalg import svd, det
from Bio.PDB import MMCIFParser


def _kabsch_Rt(Pp, Qq):
    Pc, Qc = Pp.mean(0), Qq.mean(0)
    Hh = (Pp - Pc).T @ (Qq - Qc)
    U, S, Vt = svd(Hh)
    R = Vt.T @ np.diag([1, 1, np.sign(det(Vt.T @ U.T))]) @ U.T
    return R, Qc - Pc @ R.T


def _consensus(fit, iters=12, tol=1e-4):
    ref = fit[0] - fit[0].mean(0)
    for _ in range(iters):
        al = np.empty_like(fit)
        for s in range(len(fit)):
            R, t2 = _kabsch_Rt(fit[s], ref)
            al[s] = fit[s] @ R.T + t2
        nref = al.mean(0)
        nref -= nref.mean(0)
        if np.sqrt(((nref - ref) ** 2).mean()) < tol:
            ref = nref
            break
        ref = nref
    return [_kabsch_Rt(fit[s], ref) for s in range(len(fit))]


def global_rmsf(cif_of, bait, binder_preds):
    """Per-residue all-heavy-atom RMSF over binder bait chains.

    Reference-free consensus superposition (Theseus-style), on the heavy atoms every
    binder model has. Returns (residue_ids, rmsf_values) as numpy arrays; empty with
    fewer than 2 binders.
    """
    parser = MMCIFParser(QUIET=True)
    sa_named, res_names = [], []
    for pr in binder_preds:
        ch = parser.get_structure(pr, str(cif_of[pr]))[0][bait["chain_ids"][pr]]
        d, rn = {}, {}
        for r in ch:
            if r.id[0] != " ":
                continue
            for at in r:
                if at.name.startswith("H"):
                    continue
                d[(r.id[1], at.name)] = at.coord.astype(float)
                rn.setdefault(r.id[1], set()).add(at.name)
        sa_named.append(d)
        res_names.append(rn)

    if len(sa_named) < 2:
        return np.array([], int), np.array([], float)

    cols = []
    for rid in range(1, bait["length"] + 1):
        names = set(res_names[0].get(rid, set()))
        for rn in res_names[1:]:
            names &= rn.get(rid, set())
        for n in sorted(names):
            cols.append((rid, n))

    full = np.array([[sa[c] for c in cols] for sa in sa_named])   # (n_binders, A, 3)
    col_res = np.array([rid for rid, _ in cols])
    Rts = _consensus(full)
    aligned = np.stack([full[s] @ Rts[s][0].T + Rts[s][1] for s in range(len(full))])
    peratom = np.sqrt(((aligned - aligned.mean(0)) ** 2).sum(-1).mean(0))
    g_res = np.array(sorted(set(col_res)))
    g_rmsf = np.array([np.sqrt((peratom[col_res == rid] ** 2).mean()) for rid in g_res])
    return g_res, g_rmsf
