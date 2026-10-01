import numpy as np
import pytest

from alphasift import analyze, call_patches, ptm_func, confidence_signal
from alphasift.readers import discover, load_prediction, load_sequences


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
    from alphasift.figure import figure
    out = tmp_path / "m.png"
    figure(analyze(screen, P=200), out)
    assert out.exists() and out.stat().st_size > 5000


# ------------------------------------------------- the AlphaBridge source (AF3)
def test_alphabridge_source_reads_the_sweep(screen):
    """Each bait residue takes the highest cut-off at which it is still in an
    interface; residues never in one stay at 0, the same 'absent' encoding the
    PAE channel uses."""
    from conftest import add_alphabridge
    add_alphabridge(screen)
    res = analyze(screen, P=200, source="alphabridge")

    assert res.source == "alphabridge"
    assert res.preds == ["p1", "p2", "p3"]          # only these have an AlphaBridge run
    assert res.skipped and all("AlphaBridge" in why for why in res.skipped.values())

    row = dict(zip(res.preds, res.M_signal))
    assert row["p1"][9:14] == pytest.approx(0.9)    # residues 10-14 at cut-off 0.9
    assert row["p2"][9:14] == pytest.approx(0.7)
    assert row["p1"][24] == 0.0                     # residue 25: never in an interface
    assert np.isnan(res.M_pae).all()                # this channel has no PAE of its own


def test_alphabridge_source_refuses_when_there_is_no_alphabridge_run(screen):
    """AlphaFold2 screens have no such output; the failure must say so, not
    silently produce an empty matrix."""
    with pytest.raises(FileNotFoundError, match="AlphaBridge"):
        analyze(screen, P=200, source="alphabridge")


def test_unknown_source_is_rejected(screen):
    with pytest.raises(ValueError, match="source must be"):
        analyze(screen, P=200, source="distance")


# ------------------------------------- recurrence: the width band (2026-10-01)
def test_wide_cluster_is_compared_with_wide_null_clusters():
    """No convergence by construction (each partner at an independent random
    place). The old p-value fell back to null clusters of EVERY width when fewer
    than 50 had a similar width, so this 24-residue cluster (19-42) looked
    exceptional and was called recurrent; widening the band step by step
    compares it with the widest null clusters instead, and nothing is called."""
    rng = np.random.default_rng(34)
    M = np.zeros((30, 200))
    for i in range(30):
        a = int(rng.integers(0, 188)); w = int(rng.integers(4, 12))
        M[i, a:a + w] = rng.uniform(0.5, 0.95)
    assert not [s for s in call_patches(M, P=300, rng_seed=0) if s["recurrent"]]


# --------------------------------- strength at partner level (2026-09-30)
def test_strength_counts_the_partner_not_each_piece():
    """One partner with 2 + 2 confident residues, 5 apart: it holds 4 on the
    bait, so both 2-residue pieces are strength regions. Per piece ("group",
    the earlier rule) neither has 4, so nothing is called."""
    M = np.zeros((3, 60))
    M[0, [10, 11, 18, 19]] = 0.95
    part = [(s["res_start"], s["res_end"]) for s in call_patches(M, P=100)
            if s["strength_significant"]]
    assert part == [(11, 12), (19, 20)]
    assert not [s for s in call_patches(M, P=100, single_scope="group")
                if s["strength_significant"]]


def test_partner_strength_allows_single_residue_pieces():
    M = np.zeros((2, 40))
    M[0, [5, 6, 7, 30]] = 0.9          # 4 confident residues; 31 stands alone
    sites = call_patches(M, P=100)
    assert (31, 31) in [(s["res_start"], s["res_end"]) for s in sites if s["strength_significant"]]


# ----------------------------------- 3-D merge: partner rule + attach (2026-09-30)
def _merge_case(sites, binds, positions, n_models=4):
    """sites: {name: (a, b)}; binds: {partner: {residue: signal}};
    positions: {residue: (x, y, z)} -- the same bait conformation in every model."""
    import pandas as pd
    from alphasift.interfaces import merge_interfaces
    tbl = pd.DataFrame([dict(res_start=a, res_end=b, patch=n, channel="S", significant=True)
                        for n, (a, b) in sites.items()])
    atoms = {r: np.array([xyz], float) for r, xyz in positions.items()}
    preds = [f"m{k}" for k in range(n_models)]
    ifaces, _ = merge_interfaces(tbl, {pr: atoms for pr in preds}, preds,
                                 merge_method="partner", partner_sig=binds)
    return sorted(sorted(map(tuple, r)) for r in ifaces["segments"])


def _line(a, b, x0, y=0.0):
    return {r: (x0 + 3.0 * (r - a), y, 0.0) for r in range(a, b + 1)}


def test_partner_rule_links_touching_sites_that_share_partners():
    """A and B touch and the same four partners bind both: one interface.
    C is bound by the same partners but far away: its own interface."""
    pos = {**_line(10, 13, 0), **_line(50, 53, 0, y=5.0), **_line(80, 83, 0, y=100.0)}
    binds = {f"k{i}": {r: 0.9 for r in (10, 11, 50, 51, 80, 81)} for i in range(4)}
    got = _merge_case({"A": (10, 13), "B": (50, 53), "C": (80, 83)}, binds, pos)
    assert got == [[(10, 13), (50, 53)], [(80, 83)]]


def test_partner_rule_is_two_sided():
    """B's two partners both bind A (2/2 from B's side), but they are 2 of A's
    12 (0.17 < 0.2 from A's side): touching is not enough, they stay apart --
    the TRAF2 coiled coil vs C-terminal domain case."""
    pos = {**_line(10, 13, 0), **_line(50, 53, 0, y=5.0)}
    binds = {f"k{i}": {10: 0.9} for i in range(12)}
    for k in ("k0", "k1"):
        binds[k][50] = 0.9
    got = _merge_case({"A": (10, 13), "B": (50, 53)}, binds, pos)
    assert got == [[(10, 13)], [(50, 53)]]
    binds = {k: v for k, v in binds.items() if k not in ("k10", "k11")}   # 2/10 = 0.2
    assert _merge_case({"A": (10, 13), "B": (50, 53)}, binds, pos) == [[(10, 13), (50, 53)]]


def test_short_pieces_attach_group_or_stand_alone():
    """D (2 residues) touches A and its own confident partner binds A: attached.
    E has no partner at 0.5 but sits next to A in every model: attached by
    closeness. F and G are short, far from A, touch each other and share their
    own partner: one interface of their own. H is short, alone, and shares
    nothing: its own interface."""
    pos = {**_line(10, 13, 0), **_line(30, 31, 0, y=5.0), **_line(40, 40, 6, y=5.0),
           **_line(60, 61, 0, y=100.0), **_line(64, 65, 0, y=104.0),
           **_line(90, 91, 0, y=-200.0)}
    binds = {f"k{i}": {10: 0.9, 11: 0.9} for i in range(4)}
    binds["k0"][30] = 0.9                     # D's own partner, also binds A
    binds["k5"] = {40: 0.4}                   # E: only a weak partner
    binds["k6"] = {60: 0.9, 64: 0.9}          # F and G share their own partner
    binds["k7"] = {90: 0.9}                   # H
    got = _merge_case({"A": (10, 13), "D": (30, 31), "E": (40, 40), "F": (60, 61),
                       "G": (64, 65), "H": (90, 91)}, binds, pos)
    assert got == [[(10, 13), (30, 31), (40, 40)], [(60, 61), (64, 65)], [(90, 91)]]


def test_average_linkage_is_still_available(screen):
    res = analyze(screen, P=200, merge_method="average")
    assert all(s.get("interface") for s in res.sites)


def test_alphabridge_source_runs_the_partner_merge(screen):
    """The AlphaBridge source keeps working end to end: the merge reads which
    partners bind a site from the same matrix rows (sweep cut-offs there)."""
    from conftest import add_alphabridge
    add_alphabridge(screen)
    res = analyze(screen, P=200, source="alphabridge")
    assert res.params["merge_method"] == "partner"
    assert res.sites and all(s.get("interface") for s in res.sites)
    assert len(res.interfaces) >= 1
