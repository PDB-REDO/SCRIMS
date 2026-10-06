"""
Partner annotation: for each 3-D interface, every partner that binds it
confidently, with the partner residues that bind it, tiered by confidence.

Evidence = the confidence matrix's own pairs (`alphasift.matrix.contact_pairs`:
Cb-Cb <= 8 A, direction-symmetrised PAE) seen from the partner chain, with
s = 1/(1+(PAE/d0)^2) as in the matrix. Tiers per partner residue, by its best
pair: H s >= 0.8 (PAE <= 5 A), M 0.5 <= s < 0.8, L s < 0.5. Measured on the
viz_int gold set (48 solved complexes): precision H 72 %, M 37 %, L 14 % -- the
tier is how far to trust a residue.

`short_name`, `_tier`, `_partner_annotation` and `build_partner_annotations` are
extracted verbatim from viz_int.interfaces (2026-09-29).
"""

from __future__ import annotations

import re

import numpy as np

from alphasift.structures import aa1


def short_name(n):
    """Display name for a prediction (the interactor part)."""
    return n.split('__vs__')[-1] if '__vs__' in n else re.sub(r'^fold_\d{4}_\d{2}_\d{2}_[^_]+_', '', n)


def _tier(s, tau_in, tau_high):
    return 'H' if s >= tau_high else ('M' if s >= tau_in else 'L')


def _partner_annotation(iface, iface_res, pairs, int_seq, d0_signal, tau_in, tau_high):
    pset = np.array(sorted(iface_res[iface]))
    out = []
    for pr, pp in pairs.items():
        if pp is None:
            continue
        b, j, pae = pp
        m = np.isin(b, pset)
        if not m.any():
            continue
        b, j, pae = b[m], j[m], pae[m]
        s = 1.0 / (1.0 + (pae / d0_signal) ** 2)
        conf = s >= tau_in
        if not conf.any():              # who is listed: >= 1 confident pair into the interface
            continue
        best = {}                       # each partner residue's best pair
        for jj, ss in zip(j.tolist(), s.tolist()):
            if ss > best.get(jj, -1.0):
                best[jj] = ss
        tiers = {r: _tier(v, tau_in, tau_high) for r, v in best.items()}
        runs = []                       # exact runs of consecutive residues sharing a tier
        for r in sorted(tiers):
            if runs and r == runs[-1][1] + 1 and tiers[r] == runs[-1][2]:
                runs[-1][1] = r
            else:
                runs.append([r, r, tiers[r]])
        sm = int_seq.get(pr, {})
        ranges = [[int(a), int(e), t, ''.join(aa1(sm[x]) if x in sm else 'X' for x in range(a, e + 1))]
                  for a, e, t in runs]
        tl = list(tiers.values())
        out.append({'prediction': pr, 'interactor': short_name(pr), 'ranges': ranges,
                    'n_high': tl.count('H'), 'n_medium': tl.count('M'), 'n_low': tl.count('L'),
                    'bait_residues': sorted({int(x) for x in b[conf]}),
                    'n_pairs': int(conf.sum()), 'best_pae': float(pae.min()),
                    'median_pae': float(np.median(pae[conf]))})
    out.sort(key=lambda d: (-d['n_high'], -d['n_medium'], d['best_pae'], d['prediction']))
    return out


def build_partner_annotations(interfaces, iface_res, pairs, int_seq,
                              d0_signal=10.0, tau_in=0.5, tau_high=0.8):
    """{interface -> [per-partner annotation] | None}: every partner that binds the
    interface confidently, with the partner residues that bind it.

      * a partner is LISTED iff >= 1 of its pairs into the interface has s >= tau_in
      * its FOOTPRINT is every partner residue with a pair into the interface, each
        tiered by its best pair: 'H' s >= tau_high, 'M' tau_in <= s < tau_high,
        'L' s < tau_in. `ranges` = exact runs of consecutive same-tier residues,
        [start, end, tier, one-letter sequence]; nothing is merged across gaps.
      * `bait_residues`, `n_pairs`, `median_pae` count only pairs at s >= tau_in.

    `pairs` = {prediction: (bait_res, partner_res, pae) | None}. Returns None per
    interface when NO prediction has pairwise PAE: the partner side is unknown there,
    which is not the same as empty. Predictions mapped to None in an otherwise-available
    screen are skipped.
    """
    if not pairs or all(v is None for v in pairs.values()):
        return {iface: None for iface in interfaces['interface']}
    return {iface: _partner_annotation(iface, iface_res, pairs, int_seq, d0_signal, tau_in, tau_high)
            for iface in interfaces['interface']}
