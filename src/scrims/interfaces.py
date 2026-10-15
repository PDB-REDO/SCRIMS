"""
Step 3: merge called sites into 3-D interfaces.

Two sites far apart in sequence can be one binding surface in space. For every
pair of sites, `support` is the fraction of partner models in which the two
sites' bait heavy atoms come within `d_merge` A of each other. With the default
`merge_method="partner"`, two sites LINK when they touch (support >=
merge_support) and enough of the same partners bind both; interfaces are the
connected components, and short or weakly bound sites left alone then group
with each other and attach to one neighbouring interface (`_partner_merge`).
Any scipy linkage name (e.g. "average", the rule before 2026-09-30) instead
clusters sites on 1 - support, cut at 1 - merge_support. Each site then
carries the interface it belongs to.

`merge_interfaces`, `_partner_merge`, `_components` and `patch_atoms` are
extracted verbatim from viz_int.interfaces / viz_int.structures (re-synced
2026-10-01, IMPROVEMENTS #32 and its addendum). What is NOT here, deliberately: the
partner-side contact, peptide, register and motif steps that followed them in
the source pipeline. Those depend on AlphaBridge link confidence and on the
motif machinery, both outside SCRIMS's scope.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster


def patch_atoms(bait_atoms, pr, a, b):
    """Stacked heavy atoms of bait residues a..b (inclusive) in prediction `pr`."""
    xs = [bait_atoms[pr][rid] for rid in range(a, b + 1)
          if rid in bait_atoms[pr] and len(bait_atoms[pr][rid])]
    return np.vstack(xs) if xs else np.empty((0, 3))


def _components(n, edges):
    par = list(range(n))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x
    for i, j in edges:
        par[find(i)] = find(j)
    return [find(i) for i in range(n)]


def _partner_merge(pl, support, partner_sig, merge_support, min_shared, share_small, share_large,
                   tau_in, tau_strong, tau_weak, min_len, prox_support=0.9):
    """Partner rule + attach step (merge_method='partner'; 2026-09-30, user-approved configuration B).

    LINK two patches when they touch (support >= merge_support) AND >= `min_shared` partners bind both
    (signal >= `tau_in` on a patch residue) AND shared / partners of the smaller >= `share_small` AND
    shared / partners of the larger >= `share_large`; interfaces = connected components (one link is
    enough). No partner relatedness is used.
    WAITING patches = left alone after linking AND short (< `min_len` residues) or weakly bound (< `min_shared`
    partners at `tau_in`). Their "own" partners: those at >= `tau_strong` on them if any, else at >= `tau_weak`.
    GROUP (2026-09-30, user: waiting pieces must be able to join each other): waiting patches join each other
    when they touch and share >= 1 own partner, or -- when one of them has NO partner at `tau_in` (no reliable
    partner list to compare) -- when they are close: support >= `prox_support`. Groups = connected components
    among waiting patches only.
    ATTACH: each group joins ONE real interface (never a waiting patch) when a member touches it and one of the
    group's own partners binds that neighbour's patch at >= `tau_in`, or a partner-less member is close to it
    (support >= `prox_support`); best neighbour = highest support; otherwise the group is its own interface.
    A group never links two interfaces.

    Why (validation 56_-64_, VALIDATION_LOG "§2 calling + §3 merge audit"): average linkage averaged support
    over ALL patch pairs, so a large spread-out interface diluted real neighbours (USP7 I9/I10 touch in 100 %
    of models, all I9 partners bind I10, averaged support 0.28); geometry-only single linkage chains whole
    proteins together (123 -> 72 interfaces); the one-sided partner rule let 5 TRAF-family partners glue
    TRAF2's coiled coil to its C-terminal domain (5/38 from the large side). Gold set: distinct sites
    wrongly joined 17/45 -> 5/45, same-site pairs kept together 32/39 -> 33/39 (with partner-level calling).
    Grouping + closeness (validation 67_-69_): USP7's GMPS pieces 800/809-810/829-832 and the touching weak
    patches 841-851/866-872 were stuck as five tiny interfaces; now one (800-872). Gold unchanged (5/45, 33/39),
    interfaces 135 -> 123 over 25 baits.
    """
    n = len(pl)

    def pset(i, tau):
        a, b = pl[i][0], pl[i][1]
        return {k for k, sig in partner_sig.items()
                if max((sig.get(r, 0.0) for r in range(a, b + 1)), default=0.0) >= tau}
    P5 = [pset(i, tau_in) for i in range(n)]
    P8 = [pset(i, tau_strong) for i in range(n)]
    P3 = [pset(i, tau_weak) for i in range(n)]
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            sh = len(P5[i] & P5[j])
            if (support[i, j] >= merge_support and sh >= min_shared
                    and sh / max(1, min(len(P5[i]), len(P5[j]))) >= share_small
                    and sh / max(1, max(len(P5[i]), len(P5[j]))) >= share_large):
                edges.append((i, j))
    lab0 = _components(n, edges)
    members = {l: lab0.count(l) for l in set(lab0)}
    size = {}
    for i, l in enumerate(lab0):
        size[l] = size.get(l, 0) + pl[i][1] - pl[i][0] + 1
    alone = [i for i in range(n) if members[lab0[i]] == 1]
    cand = {i for i in alone if size[lab0[i]] < min_len or len(P5[i]) < min_shared}
    own = {i: (P8[i] if P8[i] else P3[i]) for i in cand}
    cl = sorted(cand)
    glab = _components(len(cl), [(a, b) for a in range(len(cl)) for b in range(a + 1, len(cl))
                                 if (support[cl[a], cl[b]] >= merge_support and own[cl[a]] & own[cl[b]])
                                 or ((not P5[cl[a]] or not P5[cl[b]]) and support[cl[a], cl[b]] >= prox_support)])
    groups = {}
    for k, i in enumerate(cl):
        groups.setdefault(glab[k], []).append(i)
    lab = list(lab0)
    for G in groups.values():
        gown = set().union(*[own[i] for i in G])
        opts = [(support[i, j], j) for i in G for j in range(n) if j not in cand and lab0[j] != lab0[i]
                and ((support[i, j] >= merge_support and gown & P5[j])
                     or (not P5[i] and support[i, j] >= prox_support))]
        target = lab0[max(opts)[1]] if opts else lab0[G[0]]
        for i in G:
            lab[i] = target
    names = {}
    return [names.setdefault(l, len(names) + 1) for l in lab]


def merge_interfaces(patches, bait_atoms, binder_preds,
                     merge_3d=True, merge_method='average', merge_support=0.50,
                     d_merge=8.0, partner_sig=None, min_shared=2, share_small=0.7,
                     share_large=0.2, tau_in=0.5, tau_strong=0.8, tau_weak=0.3, min_len=4,
                     prox_support=0.9):
    """Merge called patches into 3-D interfaces (nb06 §3).

    For every pair of patches, `support` = fraction of binder folds in which the
    two patches' bait Cβ/heavy atoms come within `d_merge` Å. `merge_method`:
      'average' (and other scipy linkages; the geometry channel): patches clustered on
          `1 - support`, cut at `1 - merge_support`.
      'partner' (the confidence channel since 2026-09-30): the partner rule + attach step of
          `_partner_merge`; needs `partner_sig` = {binder: {bait residue: signal}}.
    Returns `(interfaces_df, iface_res)`.
    """
    pl = [(int(r.res_start), int(r.res_end), r.patch, r.channel, bool(r.significant))
          for _, r in patches.iterrows()]
    nP = len(pl)

    if merge_3d and nP > 1:
        cache = {(k, i): patch_atoms(bait_atoms, k, pl[i][0], pl[i][1])
                 for k in binder_preds for i in range(nP)}
        support = np.eye(nP)
        for i in range(nP):
            for j in range(i + 1, nP):
                dm = [float(cdist(cache[(k, i)], cache[(k, j)]).min())
                      for k in binder_preds if len(cache[(k, i)]) and len(cache[(k, j)])]
                if dm:
                    support[i, j] = support[j, i] = float((np.array(dm) < d_merge).mean())
        if merge_method == 'partner':
            if partner_sig is None:
                raise ValueError("merge_method='partner' needs partner_sig")
            comp = np.array(_partner_merge(pl, support, partner_sig, merge_support, min_shared,
                                           share_small, share_large, tau_in, tau_strong, tau_weak,
                                           min_len, prox_support))
        else:
            Dmat = 1.0 - support
            np.fill_diagonal(Dmat, 0.0)
            Dmat = (Dmat + Dmat.T) / 2
            Z = linkage(squareform(Dmat, checks=False), method=merge_method)
            comp = fcluster(Z, t=1.0 - merge_support, criterion='distance')
    else:
        comp = np.arange(nP)

    interfaces = []
    for c in sorted(set(comp)):
        mem = [i for i in range(nP) if comp[i] == c]
        segs = sorted((pl[i][0], pl[i][1]) for i in mem)
        resset = sorted({r for a, b in segs for r in range(a, b + 1)})
        interfaces.append({'interface': f'I{c}', 'n_patches': len(mem), 'segments': segs,
                           'residues': resset, 'channels': '+'.join(sorted({pl[i][3] for i in mem})),
                           'significant': any(pl[i][4] for i in mem)})
    interfaces = (pd.DataFrame(interfaces)
                  .sort_values(['significant', 'n_patches'], ascending=[False, False])
                  .reset_index(drop=True))
    iface_res = dict(zip(interfaces['interface'], interfaces['residues']))
    return interfaces, iface_res
