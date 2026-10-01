"""
Row ordering for the matrix figure: weighted-Jaccard similarity between partner
profiles, average linkage, leaf order. Extracted verbatim from
viz_int.plotdata (`_sim_wjac`, `binder_clustering`); the RMSF machinery in that
module is out of scope here and was left behind, which is why this file needs
no Bio.PDB.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from scipy.spatial.distance import squareform


def _sim_wjac(X):
    n = len(X)
    Sm = np.ones((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            mn = np.minimum(X[i], X[j]).sum()
            mx = np.maximum(X[i], X[j]).sum()
            Sm[i, j] = Sm[j, i] = (mn / mx) if mx > 0 else 0.0
    return Sm


def binder_clustering(M):
    """Return (order, clusters) for the binder × residue heatmap (nb06 ov_2).

    order    : leaf order (indices into M's rows) from average-linkage clustering of
               Gaussian-smoothed (σ=5) weighted-Jaccard profiles.
    clusters : flat cluster label per binder (inconsistency cut, maxclust=4 fallback).
    Degenerate cases (0/1 binder) return trivial order/clusters.
    """
    n = M.shape[0]
    if n == 0:
        return np.array([], int), np.array([], int)
    if n == 1:
        return np.array([0]), np.array([1])
    prof = gaussian_filter1d(M.astype(float), sigma=5, axis=1)
    Sm = _sim_wjac(prof)
    Dm = 1 - Sm
    np.fill_diagonal(Dm, 0.0)
    Dm = (Dm + Dm.T) / 2
    Z = linkage(squareform(Dm, checks=False), method="average")
    order = leaves_list(Z)
    clusters = fcluster(Z, t=1.0, depth=2, criterion="inconsistent")
    if np.unique(clusters).size < 2:
        clusters = fcluster(Z, t=4, criterion="maxclust")
    return order, clusters


