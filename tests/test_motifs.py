"""The partner side, the motif stage and the outputs drawn from them."""
import json
import shutil

import numpy as np
import pytest

from alphasift import analyze
from alphasift.matrix import bait_confidence_row, contact_pairs
from alphasift.motif_align import align_motif
from alphasift.readers import load_prediction
from alphasift.readers.mmcif import chain_identities

needs_dssp = pytest.mark.skipif(shutil.which("mkdssp") is None, reason="mkdssp (DSSP 4.x) not installed")



# ----------------------------------------------------------------- contact pairs
def test_pairs_reduce_to_the_matrix_row(screen):
    """The partner side reads the matrix's own pairs: their minimum PAE per bait
    residue IS the matrix value, so the two sides can never disagree."""
    p = load_prediction(screen, "p1")
    b, _j, pae = contact_pairs(p.chain_ids, p.res_ids, p.cb, p.pae, "A")
    rids, raw, _ = bait_confidence_row(p.chain_ids, p.res_ids, p.cb, p.pae, "A")
    for r, v in zip(rids, raw):
        got = pae[b == r]
        assert (np.isnan(v) and got.size == 0) or got.min() == v


# --------------------------------------------------------------- the run itself
@pytest.fixture(params=["dssp", "none"])
def run(screen, request):
    """One analysis with the motif stage, with and without secondary structure in the
    null (few permutations: these tests check plumbing and shapes, not calibration)."""
    if request.param == "dssp" and shutil.which("mkdssp") is None:
        pytest.skip("mkdssp (DSSP 4.x) not installed")
    return analyze(screen, P=200, motif_n_perm=50, motif_ss=request.param)


def test_partners_are_listed_with_tiers(run):
    """p1 binds bait 10-14 with its residues 1-5 at PAE 1.0 (s = 0.99): one H run."""
    site = next(i for i, v in run.partners.items() if v and 12 in run.iface_res[i])
    p1 = next(a for a in run.partners[site] if a["prediction"] == "p1")
    assert p1["ranges"] == [[1, 5, "H", "LLLLL"]]
    assert p1["n_high"] == 5 and p1["n_low"] == 0
    # p5 binds elsewhere (bait 30-34) and is not listed here
    assert "p5" not in {a["prediction"] for a in run.partners[site]}


def test_motif_engine_runs_on_every_interface(run):
    eng = run.motifs
    assert eng["engine"] == "c"
    assert set(eng["ifaces"]) == set(run.interfaces["interface"])
    site = next(i for i, v in run.partners.items() if v and 12 in run.iface_res[i])
    cols = {c["pos"]: c for c in eng["ifaces"][site]["panel1"]}
    # four partners, four different residue types -> four groups touch residue 12
    assert cols["12"]["n_groups"] == 4
    assert cols["12"]["agree"] is not None          # >= 3 groups: tested
    # every partner is its own group: no two share 40 % identity
    assert len(set(eng["groups"].values())) == len(eng["groups"])


def test_motif_stage_can_be_skipped(screen):
    res = analyze(screen, P=200, motif_engine=None)
    assert res.motifs == {}
    assert res.partners                              # the partner side still runs


def test_without_dssp_no_external_program_is_needed(screen, monkeypatch):
    """motif_ss="none": letters shuffled across the whole interface; runs with nothing on PATH."""
    import alphasift.motifs as m
    monkeypatch.setattr(m.shutil, "which", lambda name: None)
    res = analyze(screen, P=200, motif_n_perm=20, motif_ss="none")
    assert res.motifs["engine"] == "c" and res.motifs["ifaces"]
    assert res.motifs["params"]["SS"] == "none"


def test_dssp_is_the_default_and_missing_mkdssp_is_a_clear_error(screen, monkeypatch):
    from alphasift import DEFAULTS
    import alphasift.motifs as m
    assert DEFAULTS["motif_ss"] == "dssp"
    monkeypatch.setattr(m.shutil, "which", lambda name: None)
    with pytest.raises(RuntimeError, match="mkdssp.*motif-ss none"):
        analyze(screen, P=200, motif_n_perm=20)


@needs_dssp
def test_the_two_nulls_differ_only_in_how_letters_are_pooled(screen):
    """Same partners, same positions, same groups either way; only the shuffle changes."""
    a = analyze(screen, P=200, motif_n_perm=50, motif_ss="dssp").motifs
    b = analyze(screen, P=200, motif_n_perm=50, motif_ss="none").motifs
    assert a["groups"] == b["groups"] and set(a["ifaces"]) == set(b["ifaces"])
    for i in a["ifaces"]:
        pa = [(c["pos"], c["n_groups"], c["agree"]) for c in a["ifaces"][i]["panel1"]]
        pb = [(c["pos"], c["n_groups"], c["agree"]) for c in b["ifaces"][i]["panel1"]]
        assert pa == pb                                  # observed agreement does not depend on the null


# ------------------------------------------------------------- outputs
def test_document_tables_figures_and_chimerax(run, tmp_path):
    import matplotlib
    matplotlib.use("Agg")
    from alphasift import views
    from alphasift.chimerax import write_cxc
    from alphasift.export import SCHEMA, document, write_json, write_tables
    doc = document(run)
    write_json(doc, tmp_path / "a.json")
    back = json.loads((tmp_path / "a.json").read_text())
    assert back["schema"] == SCHEMA
    assert back["matrix"]["row_order"] and len(back["matrix"]["signal"]) == 6
    assert all(i["motif"] is not None for i in back["interfaces"] if i["partners"])
    write_tables(back, tmp_path)
    for f in ("partners.tsv", "bait_side.tsv", "motifs.tsv"):
        assert (tmp_path / f).exists()
    views.overview(back, tmp_path / "o.png")
    for k in range(len(back["interfaces"])):
        views.interface_card(back, k, tmp_path / f"c{k}.png")
        views.partner_strips(back, k, tmp_path / f"p{k}.png")
    assert (tmp_path / "o.png").stat().st_size > 5000
    from alphasift.report import html
    page = html(back)
    assert page.startswith("<!doctype html>") and 'id="alphasift-data"' in page
    embedded = page.split('id="alphasift-data">', 1)[1].split("</script>", 1)[0]
    assert json.loads(embedded.replace("<\\/", "</")) == back      # the page carries the whole document
    cxc = write_cxc(back, lambda pr: run.screen_dir / pr / f"{pr}_model.cif", tmp_path)
    text = cxc.read_text()
    assert "open models/" in text and "name I" in text


# ------------------------------------------------ the motif logo (website port)
def _entry(orientation="parallel"):
    return {"elements": [{"orientation": orientation, "carriers": [
        {"id": "a", "group": "a", "map": {"10": 100, "12": 103}},
        {"id": "b", "group": "b", "map": {"12": 203, "10": 200}},
        {"id": "c", "group": "c", "map": {"12": 52, "14": 70}}]}]}


def _partners():
    return [{"id": "a", "ranges": [[100, 103, "H", "LxxL"]]},
            {"id": "b", "ranges": [[200, 203, "H", "IxxL"]]},
            {"id": "c", "ranges": [[52, 52, "H", "L"], [70, 70, "M", "F"]]}]


def test_logo_anchors_on_the_position_most_carriers_hold():
    """Bait 12 is held by all three carriers (10 by two, 14 by one): it is the anchor,
    and every carrier's residue there sits at offset 0."""
    al = align_motif(_entry(), {"elements": [0], "carriers": ["a", "b", "c"]}, _partners())
    assert al["anchor"] == 12 and al["anchor_carriers"] == 3
    col0 = next(c for c in al["columns"] if c.get("offset") == 0)
    assert col0["k"] == 3 and col0["counts"] == [["L", 3]]
    colm3 = next(c for c in al["columns"] if c.get("offset") == -3)
    assert colm3["k"] == 2 and sorted(colm3["counts"]) == [["I", 1], ["L", 1]]


def test_logo_collapses_long_empty_runs():
    """c's second residue is 18 along the partner: 17 empty offsets become one gap column."""
    al = align_motif(_entry(), {"elements": [0], "carriers": ["a", "b", "c"]}, _partners())
    gaps = [c for c in al["columns"] if c["kind"] == "gap"]
    assert gaps == [{"kind": "gap", "skipped": 17}]


def test_logo_reads_antiparallel_carriers_backwards():
    al = align_motif(_entry("antiparallel"), {"elements": [0], "carriers": ["a", "b", "c"]}, _partners())
    assert al["orientation"] == "antiparallel" and all(c["dir"] == -1 for c in al["carriers"])
    assert any(c.get("offset") == 3 for c in al["columns"])     # 100 -> +3 when read C->N


# ------------------------------------------------- identities from the header
HEADER = """data_x
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'Bait protein'
2 polymer 'Prey protein'
#
loop_
_entity_poly.entity_id
_entity_poly.pdbx_seq_one_letter_code
1 MAAAK
2 MLLLK
#
loop_
_struct_ref.id
_struct_ref.db_code
_struct_ref.entity_id
_struct_ref.pdbx_db_accession
1 BAIT_HUMAN 1 P11111
2 PREY_HUMAN 2 Q22222
#
loop_
_struct_asym.id
_struct_asym.entity_id
A 2
B 1
#
"""


def test_identities_follow_sequence_not_the_header_chain_link():
    """_struct_asym says chain A is entity 2, but chain A's sequence is entity 1's:
    the sequence wins (the header link is wrong in ~37 % of raw NVIDIA files)."""
    got = chain_identities(HEADER, {"A": "MAAAK", "B": "MLLLK"})
    assert got["A"]["uniprot"] == "P11111" and got["B"]["uniprot"] == "Q22222"
    assert chain_identities(HEADER, {"A": "MXXXK"})["A"] is None
