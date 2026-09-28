import numpy as np
import pytest

from confine import analyze, call_patches, ptm_func, confidence_signal
from confine.readers import discover, load_prediction, load_sequences


# --------------------------------------------------------------- the transform
def test_d0_is_the_half_point():
    """0.5 must mean exactly PAE = d0. The figure's colourbar says so."""
    assert ptm_func(10.0, 10.0) == pytest.approx(0.5)
    assert ptm_func(5.0, 10.0) == pytest.approx(0.8)


def test_signal_is_zero_where_there_is_no_contact():
    """0 must mean absence, not weakness -- the statistics assume it."""
    sig = confidence_signal(np.array([1.0, np.nan, 30.0]), 10.0)
    assert sig[1] == 0.0
    assert sig[0] > sig[2] > 0.0


# ------------------------------------------------------------- the statistics
def test_planted_convergent_site_is_recurrent():
    M = np.zeros((20, 100))
    for i in range(20):
        M[i, 40:48] = 0.9
        M[i, (i * 3) % 90:(i * 3) % 90 + 5] = 0.6
    sites = call_patches(M, P=500, rng_seed=0)
    rec = [s for s in sites if s["recurrent"]]
    assert rec, "a site engaged by all 20 partners must be recurrent"
    assert any(s["res_start"] <= 44 <= s["res_end"] for s in rec)


def test_no_convergence_gives_no_recurrent_call():
    """Each partner binds its own place: nothing converges, nothing is recurrent."""
    rng = np.random.default_rng(0)
    M = np.zeros((20, 200))
    for i in range(20):
        a = rng.integers(0, 190)
        M[i, a:a + 8] = 0.9
    sites = call_patches(M, P=500, rng_seed=0)
    assert not [s for s in sites if s["recurrent"]]


def test_strength_needs_one_partner_to_hold_the_bar_alone():
    """Three partners covering 2 residues each is not a strength region; one
    partner covering 4 is."""
    spread = np.zeros((3, 60))
    spread[0, 10:12] = 0.95; spread[1, 12:14] = 0.95; spread[2, 14:16] = 0.95
    assert not [s for s in call_patches(spread, P=100) if s["strength_significant"]]

    alone = np.zeros((3, 60))
    alone[0, 10:14] = 0.95
    assert [s for s in call_patches(alone, P=100) if s["strength_significant"]]


# ------------------------------------------------------------------- the reader
def test_discover_and_read(screen):
    names = discover(screen)
    assert names == ["p1", "p2", "p3", "p4", "p5", "p6"]
    p = load_prediction(screen, "p1")
    assert p.format == "af3"
    assert p.pae.shape == (46, 46)
    assert set(load_sequences(screen, "p1")) == {"A", "B"}


def test_bait_is_found_by_sequence_without_alphabridge(screen):
    res = analyze(screen, P=200)
    assert res.bait["length"] == 40
    assert res.bait["n_predictions"] == 6
    assert set(res.bait["chain_ids"].values()) == {"A"}


# -------------------------------------------------------------- end to end
def test_end_to_end(screen):
    res = analyze(screen, P=500)
    assert res.M_signal.shape == (6, 40)
    assert not res.skipped

    # residue 12 is engaged by four partners at PAE 1.0-2.5 -> high signal
    col12 = res.M_signal[:, 11]
    assert (col12 > 0.8).sum() >= 3
    # residue 25 is engaged by nobody
    assert res.M_signal[:, 24].sum() == 0
    assert np.isnan(res.M_pae[:, 24]).all()

    sig = res.significant()
    assert sig, "the planted site should be called"
    assert any(s["res_start"] <= 12 <= s["res_end"] for s in sig)
    # every called site is placed in an interface
    assert all(s.get("interface") for s in res.sites)
    assert len(res.interfaces) >= 1


def test_single_pair_reports_no_qvalue(screen):
    """One prediction cannot support a convergence claim."""
    import shutil
    for name in ("p2", "p3", "p4", "p5", "p6"):
        shutil.rmtree(screen / name)
    res = analyze(screen, P=200)
    assert res.single_pair
    assert all(s["q_recur"] is None and not s["recurrent"] for s in res.sites)


def test_figure_renders(screen, tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    from confine.figure import figure
    out = tmp_path / "m.png"
    figure(analyze(screen, P=200), out)
    assert out.exists() and out.stat().st_size > 5000
