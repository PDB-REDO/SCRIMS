"""
Sequence-logo maths shared by the static figures (scrims.views) and mirrored in the
HTML report's script (scrims.report, LOGO_JS) -- a port of the companion website's
src/lib/logo.ts (2026-10-07), so the three draw the same logo.

CLASSIC LOGO (WebLogo-style, Schneider & Stephens 1990): column height = information
content R = log2(20) - (H + e(n)), floored at 0, with H = -sum f log2 f; letter height
= f x R. n is the number of INDEPENDENT partner groups (partners >= motif_id_group
identity count once) and f are the group-weighted frequencies.

e(n) is the EXACT small-sample correction: e(n) = log2(20) - E[H_n], with E[H_n] the
expected entropy of n letters drawn uniformly from 20. By linearity of expectation
over the 20 letters (each count c ~ Binomial(n, 1/20)):
    E[H_n] = -20 * sum_{c=1..n} P(c) (c/n) log2(c/n).
e(1) = 4.32 (one sequence carries no information), e(3) = 2.84, and e(n) tends to the
large-n approximation 19 / (2 ln2 n) -- which must NOT be used for small n (it is 4.57
bits at n = 3, above the 4.32 maximum, and empties every small logo).

STACKS: letters too small to draw (< MIN_GLYPH_PX) are not skipped while still taking
their height (that left an empty band under the visible letters). Visible letters sit
contiguously on the baseline, largest on top, and every too-small share is pooled into
ONE thin grey "other" band at the top; the column keeps its exact height.
"""

from __future__ import annotations

import math

LOG2_20 = math.log2(20)
MIN_GLYPH_PX = 3.0
_E_CACHE: dict[int, float] = {}


def small_sample_correction(n):
    """Exact e(n) in bits for n independent sequences (groups)."""
    if not n or n < 1:
        return LOG2_20
    k = int(round(n))
    if k in _E_CACHE:
        return _E_CACHE[k]
    p = 1 / 20
    pmf = (1 - p) ** k                       # P(c = 0)
    eh = 0.0
    for c in range(1, k + 1):
        pmf *= (k - c + 1) / c * (p / (1 - p))
        f = c / k
        eh -= 20 * pmf * f * math.log2(f)
    e = max(0.0, LOG2_20 - eh)
    _E_CACHE[k] = e
    return e


def info_content(freqs, n):
    """Information content (bits) of a frequency vector from n independent sequences, floored at 0."""
    if not n or n <= 0:
        return 0.0
    tot = sum(freqs)
    if tot <= 0:
        return 0.0
    H = -sum((f / tot) * math.log2(f / tot) for f in freqs if f > 0)
    return max(0.0, LOG2_20 - (H + small_sample_correction(n)))


def group_weighted(items):
    """items = [(letter, group id)]. Every group weighs 1 in total, split over its members
    present. Returns ([(letter, frequency)] descending, summing to 1; n = groups present)."""
    members = {}
    for _, g in items:
        members[g] = members.get(g, 0) + 1
    n = len(members)
    w = {}
    for l, g in items:
        w[l] = w.get(l, 0.0) + 1 / members[g] / n
    return sorted(w.items(), key=lambda kv: -kv[1]), n


def bits_axis_max(ics):
    """Classic-view y axis, normalised PER LOGO: the tallest column rounded up to the next
    0.5 bit (minimum 0.5), capped at log2(20) rounded up (4.5)."""
    m = max([0.0] + list(ics))
    return min(max(0.5, math.ceil(m * 2 - 1e-9) / 2), math.ceil(LOG2_20 * 2) / 2)


def build_stack(letters, min_px=MIN_GLYPH_PX):
    """letters = [(letter, height px)] whose heights sum to the column height.
    Returns (visible [(letter, h)] bottom -> top, smallest first; other_h px; other letters)."""
    vis = sorted([(l, h) for l, h in letters if h >= min_px], key=lambda x: x[1])
    small = [(l, h) for l, h in letters if 0 < h < min_px]
    return vis, sum(h for _, h in small), [l for l, _ in small]
