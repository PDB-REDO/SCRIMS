"""
The visualisation matrix: one row per prediction, one column per bait residue.

Per bait residue i, geometry only SELECTS candidate pairs (Cb-Cb <= d0_contact,
default 8 A); the value is then read from PAE alone as the direction-
symmetrised minimum PAE over the selected partner residues. That raw min-PAE
(Angstrom, NaN where this prediction has no contact) is what is displayed; for
statistics it passes through the pTM-style pair term
ptm_func(x, d0) = 1/(1+(x/d0)^2) with ONE fixed d0 = 10 A for every screen and
partner (the value pDockQ2 fits; Zhu et al., Bioinformatics 2023). The
transform flips direction (lower PAE -> higher signal) and bounds to (0, 1], so
0 means "no contact" and alphasift.stats runs on it unchanged.

`ptm_func`, `bait_confidence_row` and `confidence_signal` are extracted verbatim
from viz_int.confidence, minus the `calc_d0` length-dependent d0 helper: its only
consumer computed a per-pair `d0chn` that the loader discarded, since d0 has been
a fixed 10 A for every screen and partner since 2026-09-21. Dropping it also
keeps this package clear of any code resembling AlphaBridge's own implementation
of that published TM-score term.
"""

from __future__ import annotations

import numpy as np


# --------------------------------------------------------------- ipSAE-style transform
def ptm_func(x, d0):
    """Predicted-TM-score-style transform: bounded to (0,1], 1 at x=0, -> 0
    as x -> inf. Monotonically decreasing in x (PAE), so lower PAE -> higher
    confidence — the opposite direction of raw PAE, with no arbitrary cap."""
    return 1.0 / (1.0 + (np.asarray(x, dtype=float) / d0) ** 2)


# --------------------------------------------------------------- per-prediction reduction
def bait_confidence_row(chain_ids, res_ids, cb, pae, bait_chain, d0_contact=8.0):
    """One prediction's PAE-derived confidence, per bait residue.

    chain_ids, res_ids, cb, pae : aligned per-token arrays (T,)/(T,)/(T,3)/(T,T),
        as returned by contact_proxy.parse_cif_cb + load_confidences (AF3) or
        colabfold_adapter.load_entry (ColabFold).
    bait_chain : the bait's chain id in this prediction.

    Returns (bait_res_ids, raw_conf, d0chn, interactor_chain):
      bait_res_ids  (n_bait,) int, this prediction's bait residue numbers.
      raw_conf      (n_bait,) float, direction-symmetrised min PAE (Angstrom)
                    to the nearest (<=d0_contact A Cb-Cb) interactor residue;
                    NaN where the bait residue has no such contact.
      d0chn         ipSAE-style d0 from len(bait)+len(interactor) — NOT a
                    fixed constant, recomputed per prediction (chain lengths
                    vary with the binder).
      interactor_chain  the single non-bait chain id.
    """
    chain_ids = np.asarray(chain_ids)
    res_ids = np.asarray(res_ids)
    is_bait = chain_ids == bait_chain
    other = sorted(set(chain_ids[~is_bait].tolist()))
    if len(other) != 1:
        raise ValueError(f"expected exactly one non-bait chain, got {other}")
    interactor_chain = other[0]
    is_int = chain_ids == interactor_chain

    idx_b, idx_i = np.where(is_bait)[0], np.where(is_int)[0]
    cb_b, cb_i = cb[idx_b], cb[idx_i]

    diff = cb_b[:, None, :] - cb_i[None, :, :]
    dist = np.sqrt((diff ** 2).sum(-1))
    contact = dist <= d0_contact

    pae_bi = pae[np.ix_(idx_b, idx_i)]
    pae_ib = pae[np.ix_(idx_i, idx_b)].T
    pae_sym = np.minimum(pae_bi, pae_ib)

    raw_conf = np.full(len(idx_b), np.nan)
    has_contact = contact.any(axis=1)
    masked = np.where(contact, pae_sym, np.inf)
    raw_conf[has_contact] = masked[has_contact].min(axis=1)

    return res_ids[idx_b], raw_conf, interactor_chain


def confidence_signal(raw_conf, d0):
    """raw_conf (Angstrom, NaN=no contact) -> call_patches-ready signal
    (0=no contact, higher=more confident, bounded to (0,1])."""
    sig = np.zeros(len(raw_conf), dtype=float)
    ok = np.isfinite(raw_conf)
    sig[ok] = ptm_func(raw_conf[ok], d0)
    return sig

