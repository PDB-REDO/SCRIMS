"""Integration test on real AlphaFold3 screens: ezrin, moesin and radixin (data/, not in git).

Each screen is run through the full command line three ways:
  1. --source pae on plain AlphaFold3 output (a symlink copy without AlphaBridge files),
  2. --source pae with AlphaBridge files present (bait detection then reads them),
  3. --source alphabridge.
Checked: every prediction is in the matrix and none skipped; the bait-vs-bait homodimer has
a row; at least one significant site; every output file is written; the motif engine ran on
every interface with partners; 1 and 2 are identical; the AlphaBridge matrix has no PAE;
every significant residue the AlphaBridge source calls is also called by the PAE source
(true on all three screens, 2026-10-07).

Slow (~5 min, DSSP null): run with `pytest --run-slow`. Skipped when data/ is absent, and
the AlphaBridge cases when a screen has no AlphaBridge runs (AlphaBridge's
define_interfaces.py over every prediction folder creates them).
"""
import json
from pathlib import Path

import pytest

from scrims.cli import main

DATA = Path(__file__).resolve().parents[1] / "data"
SCREENS = [("EZR_screen", 61, "p15311"), ("MSN_screen", 41, "p26038"), ("RDX_screen", 24, "p35241")]
OUTPUTS = ("scrims.json", "report.html", "overview.png", "matrix.png", "sites.tsv", "interfaces.tsv",
           "partners.tsv", "bait_side.tsv", "motifs.tsv", "matrix.npz", "summary.json", "scrims.cxc")

pytestmark = pytest.mark.slow


def _screen(name):
    d = DATA / name
    if not d.is_dir():
        pytest.skip(f"{d} not present (data/ is not in git)")
    return d


def _af3_only(src, dst):
    """The screen without AlphaBridge files: per prediction, symlinks to its top-level AF3 files."""
    for pred in sorted(p for p in src.iterdir() if p.is_dir()):
        (dst / pred.name).mkdir(parents=True)
        for f in pred.glob(f"{pred.name}_*"):
            if f.is_file():
                (dst / pred.name / f.name).symlink_to(f)
    return dst


def _run(screen, out, *extra):
    assert main(["run", str(screen), "-o", str(out), *extra]) == 0
    return json.loads((out / "scrims.json").read_text()), json.loads((out / "summary.json").read_text())


def _check_run(out, doc, summary, n, bait_acc):
    assert summary["n_predictions_in_matrix"] == n and not summary["skipped"]
    assert [p["id"] for p in doc["predictions"] if p["id"].endswith(f"__vs__{bait_acc}")], "homodimer row missing"
    assert summary["n_significant"] > 0
    for f in OUTPUTS:
        assert (out / f).stat().st_size > 0, f
    assert len(list((out / "interfaces").glob("*.png"))) == 2 * len(doc["interfaces"])
    assert all(i["motif"] is not None for i in doc["interfaces"] if i["partners"])
    page = (out / "report.html").read_text()
    assert 'id="scrims-data"' in page


def _sig_residues(doc):
    return {r for s in doc["sites"] if s["significant"] for r in range(s["res_start"], s["res_end"] + 1)}


@pytest.mark.parametrize("name,n,bait_acc", SCREENS)
def test_screen_with_both_sources(name, n, bait_acc, tmp_path):
    screen = _screen(name)
    plain = _af3_only(screen, tmp_path / "af3only" / name)          # same folder name: same motif seeds
    d1, s1 = _run(plain, tmp_path / "pae_af3only")
    _check_run(tmp_path / "pae_af3only", d1, s1, n, bait_acc)
    d2, s2 = _run(screen, tmp_path / "pae")
    _check_run(tmp_path / "pae", d2, s2, n, bait_acc)
    for k in ("predictions", "matrix", "tracks", "sites", "interfaces"):
        assert d1[k] == d2[k], f"PAE results depend on AlphaBridge files being present ({k})"
    if not all((p / "AlphaBridge" / "alphabridge_data.json").exists() for p in screen.iterdir() if p.is_dir()):
        pytest.skip(f"{name}: not every prediction has an AlphaBridge run")
    d3, s3 = _run(screen, tmp_path / "alphabridge", "--source", "alphabridge")
    _check_run(tmp_path / "alphabridge", d3, s3, n, bait_acc)
    assert [p["id"] for p in d3["predictions"]] == [p["id"] for p in d2["predictions"]]
    assert all(v is None for row in d3["matrix"]["pae"] for v in row)
    assert _sig_residues(d3) <= _sig_residues(d2)
