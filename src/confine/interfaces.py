"""
Step 3: merge called sites into 3-D interfaces.

Two sites far apart in sequence can be one binding surface in space. For every
pair of sites, `support` is the fraction of partner models in which the two
sites' bait heavy atoms come within `d_merge` A of each other; sites are
clustered on 1 - support and cut at 1 - merge_support. Each site then carries
the interface it belongs to.

`merge_interfaces` and `patch_atoms` are extracted verbatim from
viz_int.interfaces / viz_int.structures. What is NOT here, deliberately: the
partner-side contact, peptide, register and motif steps that followed them in
the source pipeline. Those depend on AlphaBridge link confidence and on the
motif machinery, both outside CONFINE's scope.
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
def merge_interfaces(patches, bait_atoms, binder_preds,
                     merge_3d=True, merge_method='average', merge_support=0.50,
                     d_merge=8.0):
    """Merge called patches into 3-D interfaces (nb06 §3).

    For every pair of patches, `support` = fraction of binder folds in which the
    two patches' bait Cβ/heavy atoms come within `d_merge` Å; patches are then
    clustered on `1 - support` (`merge_method` linkage, cut at `1 - merge_support`).
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


