"""The sequence-logo rules shared by the HTML report and the figures (website round 3,
2026-10-07): exact small-sample correction, group weighting, stacking without gaps."""
import itertools
import json
import math
import shutil
import subprocess

import pytest

from scrims.logo import (LOG2_20, bits_axis_max, build_stack, group_weighted, info_content,
                            small_sample_correction)

needs_node = pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")


# ------------------------------------------------- e(n), the exact correction
def test_small_sample_correction_known_values():
    """e(1) = log2 20: one sequence carries no information; e(3) and e(6) as on the website."""
    assert small_sample_correction(1) == pytest.approx(LOG2_20, abs=1e-12)
    assert round(small_sample_correction(1), 4) == 4.3219
    assert round(small_sample_correction(3), 4) == 2.8359
    assert round(small_sample_correction(6), 4) == 1.9769
    # a conserved 3-group column: 4.32 - 2.84 = 1.49 bits (the large-n formula would give 0)
    assert info_content([1.0], 3) == pytest.approx(1.486, abs=1e-3)
    assert info_content([1.0], 1) == 0.0


def _brute_e(n):
    """log2 20 - expected entropy of n letters drawn uniformly from 20, by summing over every
    composition (count vector) with its multinomial probability."""
    eh = 0.0
    for comp in itertools.combinations_with_replacement(range(20), n):
        counts = {}
        for a in comp:
            counts[a] = counts.get(a, 0) + 1
        prob = math.factorial(n) / math.prod(math.factorial(c) for c in counts.values()) / 20 ** n
        eh += prob * -sum(c / n * math.log2(c / n) for c in counts.values())
    return LOG2_20 - eh


@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_small_sample_correction_matches_enumeration(n):
    assert small_sample_correction(n) == pytest.approx(_brute_e(n), abs=1e-9)


def test_small_sample_correction_tends_to_the_large_n_formula():
    n = 400
    assert small_sample_correction(n) == pytest.approx(19 / (2 * math.log(2) * n), rel=0.02)


def test_bits_axis_is_scaled_to_the_tallest_column():
    assert bits_axis_max([1.73, 0.2]) == 2.0       # VAPA I1 M1 on the website: axis 0-2
    assert bits_axis_max([1.49]) == 1.5
    assert bits_axis_max([0.0]) == 0.5
    assert bits_axis_max([4.32]) == 4.5


# ------------------------------------------------- group weighting
def test_group_weighted_counts_each_group_once():
    """Two paralogs (one group) carrying L and I, one unrelated partner carrying F:
    the group weighs 1 in total -> L 0.25, I 0.25, F 0.5; n = 2 groups."""
    freq, n = group_weighted([("L", "g1"), ("I", "g1"), ("F", "g2")])
    assert n == 2 and dict(freq) == pytest.approx({"L": 0.25, "I": 0.25, "F": 0.5})
    assert sum(f for _, f in freq) == pytest.approx(1.0)


# ------------------------------------------------- stacks
def test_stack_keeps_the_column_height_and_pools_small_shares():
    letters = [("L", 40.0), ("I", 20.0), ("V", 2.0), ("A", 1.0), ("F", 17.0)]
    vis, other, small = build_stack(letters)
    assert [l for l, _ in vis] == ["F", "I", "L"]            # smallest at the bottom, largest on top
    assert other == pytest.approx(3.0) and sorted(small) == ["A", "V"]
    assert sum(h for _, h in vis) + other == pytest.approx(80.0)


def test_figure_stack_has_no_gap_at_the_baseline():
    """The drawn glyph outlines sit exactly on the baseline and on each other, the grey
    band directly on top (y grows downwards in the figure canvas)."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.patches import PathPatch, Rectangle
    from scrims.views import _Canvas
    cv = _Canvas(100, 120)
    base = 100.0
    cv.stack([("Q", 30.0), ("L", 40.0), ("G", 8.0), ("W", 1.5), ("C", 0.5)], 10, base, 16)
    glyphs = sorted((p.get_path().get_extents() for p in cv.ax.patches if isinstance(p, PathPatch)),
                    key=lambda b: -b.y1)
    assert len(glyphs) == 3
    assert glyphs[0].y1 == pytest.approx(base, abs=1e-6)            # lowest letter on the baseline
    for lo, hi in zip(glyphs, glyphs[1:]):
        assert hi.y1 == pytest.approx(lo.y0, abs=1e-6)              # contiguous, no gap
    band = [p for p in cv.ax.patches if isinstance(p, Rectangle)]
    assert len(band) == 1 and band[0].get_y() + band[0].get_height() == pytest.approx(glyphs[-1].y0, abs=1e-6)
    assert base - band[0].get_y() == pytest.approx(80.0, abs=1e-6)  # the column keeps its exact height
    import matplotlib.pyplot as plt
    plt.close(cv.fig)


# ------------------------------------------------- the report's JS twin
@needs_node
def test_report_logo_js_matches_python(tmp_path):
    from scrims.report import LOGO_JS
    script = tmp_path / "t.js"
    script.write_text(LOGO_JS + """
const out = { e: [], ic: infoContent([0.5, 0.25, 0.25], 5), gw: groupWeighted([['L', 'g1'], ['I', 'g1'], ['F', 'g2']]),
  st: buildStack([['L', 40], ['I', 20], ['V', 2], ['A', 1], ['F', 17]]), ax: bitsAxisMax([1.73, 0.2]) };
for (let n = 1; n <= 40; n++) out.e.push(smallSampleCorrection(n));
console.log(JSON.stringify(out));
""")
    got = json.loads(subprocess.run(["node", str(script)], capture_output=True, text=True, check=True).stdout)
    assert got["e"] == pytest.approx([small_sample_correction(n) for n in range(1, 41)], abs=1e-12)
    assert got["ic"] == pytest.approx(info_content([0.5, 0.25, 0.25], 5), abs=1e-12)
    assert dict(got["gw"]["freq"]) == pytest.approx({"L": 0.25, "I": 0.25, "F": 0.5}) and got["gw"]["n"] == 2
    vis, other, _ = build_stack([("L", 40), ("I", 20), ("V", 2), ("A", 1), ("F", 17)])
    assert [[v["letter"], v["h"]] for v in got["st"]["visible"]] == [list(x) for x in vis]
    assert got["st"]["other"] == pytest.approx(other) and got["ax"] == 2.0
