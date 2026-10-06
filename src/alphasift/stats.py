"""
The statistical engine: recurrence FDR + strength criterion on a
partner x bait-residue matrix.

Extracted verbatim from viz_int.patches (the `call_patches` engine and its
helpers) so that the published numbers apply to this code unchanged; re-synced
2026-10-01 with viz_int's partner-level strength calling (IMPROVEMENTS #32) and
widened width band (#33). Only two things differ from the source: the legacy
geometry loaders that lived in the same module are not here (AlphaSIFT builds its
matrix from PAE, see alphasift.matrix), and the two dead parameters
`tau_strong`/`strong_mask` are gone from the signature -- they were ignored by
the body since 2026-09-21.

Recurrence asks: do many partners converge on the same bait residues more than
a positional null allows? Each partner's whole row is shifted by a random
circular offset, clusters are formed above a percentile of that null, and each
observed cluster's mass is compared against null clusters of similar extent
(the width band widens until it holds enough of them), BH-corrected.

Strength asks a different question and is NOT a test: a partner qualifies when
it holds >= `min_single` residues at >= `tau_single` anywhere on the bait, and
every sequence piece of it (gaps of up to `single_gap`, any length) is a
strength region. It is reported as a declared criterion with no p-value,
because a value claim has no null in a matrix that only supports positional
nulls.
"""

from __future__ import annotations

import numpy as np
from numpy.random import default_rng
from scipy.stats import false_discovery_control


def _extract_segments(p):
    segs, st = [], None
    for i, v in enumerate(p):
        if v > 0 and st is None:
            st = i
        elif v == 0 and st is not None:
            segs.append(p[st:i].copy()); st = None
    if st is not None:
        segs.append(p[st:].copy())
    return segs



def _footprints(M, tau, gap_max, min_res):
    """Per-partner confident footprints: for each row, its residues at >= `tau`, grouped
    so that residues separated by at most `gap_max` non-qualifying positions stay in one
    group, keeping groups of >= `min_res` residues. Returns (start, end, row) spans.

    This is the strength channel's region-forming rule (2026-09-21). It replaces a union
    across partners at a LOOSER threshold: regions used to be formed where ANY partner
    had a contact at PAE <= tau_strong_pae (10 A) and then judged on what ONE partner did
    at tau_single (5 A), so 63 of the shipped candidate regions could never pass (e.g.
    DAP3 171-174: three partners covering <= 3 residues each, best PAE 8.6 A), while a
    single partner's confident footprint whose residues were spaced more than one apart
    formed no region at all. Forming and judging now use the same evidence, so every
    strength region is significant by construction. `gap_max` is explicit because a
    partner's contacts are rarely contiguous (RHOQ 164-180: 9 partners with >= 4 confident
    residues, none of them 4 in a row).
    """
    out = []
    for i in range(M.shape[0]):
        idx = np.where(M[i] >= tau)[0]
        if len(idx) < min_res:
            continue
        st = 0
        for k in range(1, len(idx) + 1):
            if k == len(idx) or idx[k] - idx[k - 1] - 1 > gap_max:
                grp = idx[st:k]
                if len(grp) >= min_res:
                    out.append((int(grp[0]), int(grp[-1]), i))
                st = k
    return out


def _partner_footprints(M, tau, gap_max, min_res):
    """Strength footprints at PARTNER level (2026-09-30, user-approved "configuration B"): a row qualifies when it
    holds >= `min_res` residues at >= `tau` ANYWHERE on the bait, and then EVERY sequence group of it (residues at
    >= `tau`, gaps <= `gap_max`) is a footprint, whatever its length. Still sequence-level: the 3-D decision is left
    to the merge (alphasift.interfaces, partner rule + attach step).

    Why: `_footprints` asks each sequence group for `min_res` residues, which cuts one partner's footprint on a
    repeat surface into pieces and drops the short ones (ILK x LIMS4: 14 residues at s >= 0.8, only 97-101 and
    130-134 kept; 66-68, 76, 105, 109, 139 lost). Measured on the gold set (validation 57_): coverage 81 -> 86 %,
    the added residues as precise as the shipped ones (39 vs 41 %).
    """
    out = []
    for i in range(M.shape[0]):
        if int((M[i] >= tau).sum()) < min_res:
            continue
        out += [(a, b, i) for a, b, _ in _footprints(M[i:i + 1], tau, gap_max, 1)]
    return out


def _find_runs(mask, gap=1, ms=4):
    runs, st, en, gr = [], None, None, 0
    for j, x in enumerate(mask):
        if x:
            if st is None:
                st = j
            en = j; gr = 0
        elif st is not None:
            gr += 1
            if gr > gap:
                runs.append((st, en)); st = None; gr = 0
    if st is not None:
        runs.append((st, en))
    return [(a, b) for a, b in runs if b - a + 1 >= ms]


def _pmass(s, a, b):
    return float(s[a:b + 1].sum())


def call_patches(M, res_ids=None, P=5000, gap=1, min_size=4, pct=95,
                 rng_seed=0, strength_amount=None,
                 tau_single=0.8, min_single=4, single_gap=2, single_scope="partner"):
    """Call interaction sites on a partner x bait-residue matrix.

    `M` is the signal matrix (0 = no contact, higher = more confident;
    alphasift.matrix.confidence_signal produces it). Returns a list of dicts
    sorted by start: res_start, res_end, channel ('R'/'S'/'R+S'), q_recur,
    recurrent, recurrent_core, strength_A, strength_A_per_res,
    strength_significant, strength_core, n_partners_strong, significant.

    RECURRENCE (a test). Clusters of high summed signal, significance against a
    null that shifts each partner's whole row by a random circular offset --
    randomising only each partner's position relative to the others, which is
    exactly what recurrence claims. Each observed cluster is compared against
    null clusters of similar extent (+-1.5x, the band widening x1.25 per step
    until it holds >= 50 members), then BH-corrected. `pct` sets the
    cluster-forming percentile of the null.

    STRENGTH (a declared criterion, NOT a test). Per partner, residues at
    >= `tau_single` are grouped allowing gaps of up to `single_gap`.
    `single_scope` says where `min_single` applies: "partner" (default,
    2026-09-30) -- the partner needs >= `min_single` such residues anywhere on
    the bait, and then every group of it is a region, whatever its length;
    "group" (the rule of 2026-09-21) -- each group needs >= `min_single`
    residues. Overlapping groups merge across partners. Forming and judging use
    the same evidence, so every strength region is significant by construction. There is no p-value: a value claim
    ("the model is very confident here") has no null in a matrix that only
    supports positional nulls -- five designs were tested and all failed
    calibration. Report it as a criterion.

    `strength_amount` (optional, shape (L,)) is the total confident engagement
    per residue summed over all partners; it is reported as `strength_A` /
    `strength_A_per_res` and never gates significance.

    Parameters `tau_strong`, `strong_mask` and `tau_amount` exist in the source
    version and are ignored by its body (tau_amount since 2026-09-21: an A/residue
    gate next to the footprint criterion would contradict it); they are not part of
    this API.
    """
    if res_ids is None:
        res_ids = np.arange(1, M.shape[1] + 1)
    s_sum, s_max = M.sum(0), M.max(0)
    N, L = M.shape
    rng = default_rng(rng_seed)

    # Null: shift each prediction's WHOLE row by a random offset (circular), so a
    # partner's segments keep their own spacing and only its position relative to the
    # other partners is randomised -- which is the only thing recurrence claims. The
    # earlier null relocated every segment independently; measured 2026-09-21, that under-produces clusters, because a real
    # partner's segments sit close together while the null scattered them across the
    # bait. Consequence: the test was anti-conservative -- on screens built from real
    # partner rows of OTHER proteins placed at random offsets (no convergence by
    # construction) it called "recurrent" clusters in 17 of 24 screens (45 clusters);
    # with this null and the extent-matched comparison below, 2 of 24 (2 clusters),
    # while real calls went 143 -> 149 and all 16 regression sites are kept.
    # See MEMORY.md 2026-09-21 and IMPROVEMENTS.md #29. Under viz_int's deterministic
    # harness (25 baits x 5 fakes, P 5000, pct 95) the fallback-to-all-widths p-value
    # gave 12/125 (9.6 %); the width-band widening below (2026-10-01) brings it to
    # 2/125 (1.6 %) -- IMPROVEMENTS.md #33.
    ns = np.empty((P, L))
    for p in range(P):
        acc = np.zeros(L)
        for i in range(N):
            acc += np.roll(M[i], int(rng.integers(0, L)))
        ns[p] = acc
    t = float(np.percentile(ns, pct))
    # Pool keeps each null cluster's WIDTH alongside its mass: with whole rows kept
    # intact the null's clusters are wider and heavier, so comparing an observed
    # cluster against the pooled mass of clusters of every width under-calls narrow
    # real sites (measured: 9/16 regression sites kept without this, 16/16 with it).
    # Each observed cluster is therefore compared with null clusters of similar extent.
    pool = np.array([(_pmass(ns[p], a, b), b - a + 1)
                     for p in range(P) for a, b in _find_runs(ns[p] > t, gap, min_size)]
                    or [(0.0, 1)])

    # Width band (2026-10-01, user-approved "F1"): when fewer than 50 null clusters have a similar width
    # (w/1.5 .. 1.5w), WIDEN the band step by step (x1.25) until there are, instead of falling back to the
    # whole pool. The old fallback compared a wide cluster with mostly short null clusters, so every wide
    # cluster looked exceptional: fake-screen false calls 12/125 (9.6 %, all 30-68 residues wide). Widening:
    # 2/125 (1.6 %), and the 25 real screens' recurrent calls are unchanged (16 shipped runs relied on the
    # fallback; all stay significant against the widened, similar-width null). viz_int validation 66_, 70_.
    def _pval(a, b):
        m, w = _pmass(s_sum, a, b), b - a + 1
        f = 1.5
        band = (pool[:, 1] >= w / f) & (pool[:, 1] <= w * f)
        while band.sum() < 50 and f < 1e3:
            f *= 1.25
            band = (pool[:, 1] >= w / f) & (pool[:, 1] <= w * f)
        ref = pool[band, 0] if band.sum() >= 50 else pool[:, 0]
        return (int((ref >= m).sum()) + 1) / (len(ref) + 1)

    obs = _find_runs(s_sum > t, gap, min_size)
    pvals = [_pval(a, b) for a, b in obs]
    qvals = false_discovery_control(pvals, method="bh") if obs else []
    recq = {ab: q for ab, q in zip(obs, qvals)}
    # Strength regions = per-partner confident footprints (see `_footprints`), merged
    # where they overlap. (The source version's `strong_mask` detour is gone: it existed because
    # each prediction's signal used its own d0, so a scalar threshold on M meant a
    # different PAE per partner; with the fixed d0_signal (schema 0.7) M >= tau_single
    # is exactly PAE <= 5 A for every partner, so the raw-Angstrom detour was unnecessary.)
    fp = (_partner_footprints if single_scope == "partner" else _footprints)(M, tau_single, single_gap,
                                                                              min_single)
    stg = []
    for a_, b_, _ in sorted(fp):
        if stg and a_ <= stg[-1][1] + 1:
            stg[-1][1] = max(stg[-1][1], b_)
        else:
            stg.append([a_, b_])
    stg = [(int(a_), int(b_)) for a_, b_ in stg]

    items = [(a, b, "R") for a, b in obs] + [(a, b, "S") for a, b in stg]
    items.sort()
    merged = []
    for a, b, ch in items:
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b); merged[-1][2].add(ch)
        else:
            merged.append([a, b, {ch}])

    amt = None if strength_amount is None else np.asarray(strength_amount, dtype=float)

    # Strength (2026-09-21): a DECLARED confidence criterion, not a test -- one partner
    # holding >= min_single residues at >= tau_single anywhere in a strength run. It is
    # judged per RUN, not over the merged patch, so a compact confident site is not
    # diluted by whatever recurrence run it happens to overlap (measured: 6 patches
    # failed the old merged-span gate although one of their own runs passed).
    # `min_single` residues need NOT be consecutive: a partner's contacts along the bait
    # are rarely contiguous (RHOQ 164-180 has 9 partners with >= 4 residues at >= 0.8 and
    # NONE with 4 in a row), which is why an earlier "4 consecutive" form lost most of
    # the watch set. Why no p-value: with only this matrix there is no null for a VALUE
    # claim -- every version tried was either anti-conservative (63 % of null screens
    # flagged) or had zero power; see ERRORS.md 2026-09-21. Multi-partner convergence is
    # recurrence's question, so no amount/cluster term enters here (user, 2026-09-21).
    strong_runs = stg           # qualifying by construction; see `_footprints`

    out = []
    for s, e, chs in merged:
        qv = [recq[ab] for ab in obs if not (ab[1] < s or ab[0] > e)]
        q_recur = float(min(qv)) if qv else float("nan")
        recurrent = ("R" in chs) and (q_recur < 0.05)
        if amt is None:
            A = a_res = float("nan")
        else:
            A = float(amt[s:e + 1].sum()); a_res = A / (e - s + 1)
        # Evidence cores: which parts of a merged patch actually carry each call. A
        # merged patch is the union of overlapping runs, so 24 % of the residues in
        # recurrent patches used to sit outside any significant recurrence cluster.
        rec_core = [(int(res_ids[a]), int(res_ids[b])) for a, b in obs
                    if recq[(a, b)] < 0.05 and not (b < s or a > e)]
        str_core = [(int(res_ids[a]), int(res_ids[b])) for a, b in strong_runs
                    if not (b < s or a > e)]
        n_strong = len({i for a, b, i in fp if not (b < s or a > e)})
        out.append({
            "res_start": int(res_ids[s]), "res_end": int(res_ids[e]),
            "channel": "R+S" if chs == {"R", "S"} else next(iter(chs)),
            "q_recur": q_recur,
            "recurrent": recurrent,
            "recurrent_core": rec_core,
            "strength_A": A,
            "strength_A_per_res": a_res,
            "strength_significant": bool(str_core),
            "strength_core": str_core,
            "n_partners_strong": n_strong,
            "significant": recurrent or bool(str_core),
        })
    return sorted(out, key=lambda d: d["res_start"])


def significant_patches(M, res_ids=None, recurrence_only=True, **kw):
    """Return [(res_start, res_end), ...].

    recurrence_only=True (default) gives the **common patches** — residues that
    many interactors converge on (recurrence FDR<0.05), which is the project's
    definition of a hotspot. Set False for the full significant set (R ∪ S).
    """
    key = "recurrent" if recurrence_only else "significant"
    return [(p["res_start"], p["res_end"]) for p in call_patches(M, res_ids, **kw) if p[key]]
