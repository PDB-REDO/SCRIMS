"""The partner side, the motif stage and the outputs drawn from them."""
import json
import shutil

import numpy as np
import pytest

from scrims import analyze
from scrims.matrix import bait_confidence_row, contact_pairs
from scrims.motif_align import align_motif
from scrims.readers import load_prediction
from scrims.readers.mmcif import chain_identities

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


def test_bait_side_letters_grouped(run):
    """letters_grouped (the bait-side classic logo's input): per column it sums to 1, and with
    every partner its own group it equals the plain shares (vote-weighted per group)."""
    from scrims.export import document
    for r in run.motifs["ifaces"].values():
        for c in (r or {}).get("panel1", []):
            assert set(c["letters_grouped"]) == set(c["letters"])
            if c["letters_grouped"]:
                assert sum(c["letters_grouped"].values()) == pytest.approx(1.0, abs=0.01)
    doc = document(run)
    p1 = [c for i in doc["interfaces"] if i["motif"] for c in i["motif"]["panel1"]]
    assert p1 and all("letters_grouped" in c for c in p1)
    for i in doc["interfaces"]:
        for al in ((i["motif"] or {}).get("alignments") or {}).values():
            for car in al["carriers"]:
                assert set(car["map"].values()) == set(car["residues"])     # bait position -> partner residue
                assert len(car["cells"]) == (len(car["residues"]) if car["aligned"] else 0)


def test_motif_stage_can_be_skipped(screen):
    res = analyze(screen, P=200, motif_engine=None)
    assert res.motifs == {}
    assert res.partners                              # the partner side still runs


def test_without_dssp_no_external_program_is_needed(screen, monkeypatch):
    """motif_ss="none": letters shuffled across the whole interface; runs with nothing on PATH."""
    import scrims.motifs as m
    monkeypatch.setattr(m.shutil, "which", lambda name: None)
    res = analyze(screen, P=200, motif_n_perm=20, motif_ss="none")
    assert res.motifs["engine"] == "c" and res.motifs["ifaces"]
    assert res.motifs["params"]["SS"] == "none"


def test_dssp_is_the_default_and_missing_mkdssp_is_a_clear_error(screen, monkeypatch):
    from scrims import DEFAULTS
    import scrims.motifs as m
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
    from scrims import views
    from scrims.chimerax import write_cxc
    from scrims.export import SCHEMA, document, write_json, write_tables
    doc = document(run)
    write_json(doc, tmp_path / "a.json")
    back = json.loads((tmp_path / "a.json").read_text())
    assert back["schema"] == SCHEMA
    assert back["matrix"]["row_order"] and len(back["matrix"]["signal"]) == 6
    assert all(i["motif"] is not None for i in back["interfaces"] if i["partners"])
    assert all("touching" in m for i in back["interfaces"] if i["motif"] for m in i["motif"]["motifs"])
    write_tables(back, tmp_path)
    for f in ("partners.tsv", "bait_side.tsv", "motifs.tsv"):
        assert (tmp_path / f).exists()
    assert (tmp_path / "bait_side.tsv").read_text().split("\n")[0].endswith("partner_letters\tpartner_letters_grouped")
    views.overview(back, tmp_path / "o.png")
    for k in range(len(back["interfaces"])):
        views.interface_card(back, k, tmp_path / f"c{k}.png")
        views.partner_strips(back, k, tmp_path / f"p{k}.png")
    assert (tmp_path / "o.png").stat().st_size > 5000
    from scrims.report import html
    page = html(back)
    assert page.startswith("<!doctype html>") and 'id="scrims-data"' in page
    embedded = page.split('id="scrims-data">', 1)[1].split("</script>", 1)[0]
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
    """b's and c's second residue is 18 along the partner: 17 empty offsets become one gap column.
    (Two groups hold +18: since 2026-10-07 a column needs >= 2 groups.)"""
    e = _entry()
    e["elements"][0]["carriers"][1]["map"]["14"] = 221
    al = align_motif(e, {"elements": [0], "carriers": ["a", "b", "c"]}, _partners())
    gaps = [c for c in al["columns"] if c["kind"] == "gap"]
    assert gaps == [{"kind": "gap", "skipped": 17}]


def test_logo_drops_columns_held_by_one_group():
    """1b (2026-10-07): with the shared fixture only c holds +18 -- one group is 0 bits, blank space -- so it is no
    column; c keeps the residue in its cells (the per-carrier sequence still shows it)."""
    al = align_motif(_entry(), {"elements": [0], "carriers": ["a", "b", "c"]}, _partners())
    assert [c["offset"] for c in al["columns"] if c["kind"] == "col"] == [-3, -2, -1, 0]
    assert {c["id"]: c["cells"] for c in al["carriers"]}["c"] == {"0": "L", "18": "F"}


def test_logo_votes_direction_even_when_elements_agree():
    """1a (2026-10-07): VAPA I3 OSBPL7 binds antiparallel inside a 'parallel' element; read forwards, its residues
    landed 100+ offsets away. Its positions vote: c's residues fall as the bait rises, so it is read backwards."""
    entry = {"elements": [{"orientation": "parallel", "carriers": [
        {"id": "a", "group": "a", "map": {"10": 100, "12": 102, "14": 104}},
        {"id": "b", "group": "b", "map": {"10": 200, "12": 202, "14": 204}},
        {"id": "c", "group": "c", "map": {"10": 54, "12": 52, "14": 50}}]}]}
    al = align_motif(entry, {"elements": [0], "carriers": ["a", "b", "c"]}, None)
    d = {c["id"]: c["dir"] for c in al["carriers"]}
    assert d == {"a": 1, "b": 1, "c": -1}
    assert {c["id"]: sorted(int(o) for o in c["cells"]) for c in al["carriers"]}["c"] == [0, 2, 4]   # 54, 52, 50 on baits 10, 12, 14 -- as a and b


def test_logo_shows_an_antiparallel_motif_n_to_c_and_flags_it():
    """2026-10-09 (user): a motif whose partners run against the bait is MIRRORED so it reads N->C, and flagged
    (`flipped`); before, every carrier was read backwards (100 at +3)."""
    al = align_motif(_entry("antiparallel"), {"elements": [0], "carriers": ["a", "b", "c"]}, _partners())
    assert al["orientation"] == "antiparallel" and al["flipped"] is True
    assert all(c["dir"] == 1 for c in al["carriers"])
    assert any(c.get("offset") == -3 for c in al["columns"])    # 100 -> -3: read N->C
    from scrims.motif_align import mirror
    back = mirror(al)
    assert back["flipped"] is False and all(c["dir"] == -1 for c in back["carriers"])
    assert [c.get("offset") for c in back["columns"] if c["kind"] == "col"] == sorted(
        -c["offset"] for c in al["columns"] if c["kind"] == "col")


def test_logo_cells_hold_each_carriers_letter_per_offset():
    """The per-carrier table reads `cells`: a's residues 100 (L) and 103 (L) at offsets -3 and 0;
    c holds the anchor (52) and 70 (F) at +18; the absent ones are simply not there."""
    al = align_motif(_entry(), {"elements": [0], "carriers": ["a", "b", "c"]}, _partners())
    cells = {c["id"]: c["cells"] for c in al["carriers"]}
    assert cells["a"] == {"-3": "L", "0": "L"} and cells["b"] == {"-3": "I", "0": "L"}
    assert cells["c"] == {"0": "L", "18": "F"}
    assert al["carriers"][0]["map"] == {"10": 100, "12": 103}



def test_logo_register_by_consensus_not_one_position():
    """VAPA FFAT (2026-10-07): two equivalent residues (E-F-F-D-A's two F) swap on the anchor bait position in one
    carrier. b put its FIRST F on bait 52 where a and c put their second, but b's other positions agree with a and c;
    the consensus step lines b up with them instead of shifting the whole row by one."""
    entry = {"elements": [{"orientation": "parallel", "carriers": [
        {"id": "a", "group": "a", "map": {"50": 352, "52": 351, "53": 350, "56": 354, "64": 353}},
        {"id": "b", "group": "b", "map": {"52": 345, "53": 346, "56": 349, "64": 348}},
        {"id": "c", "group": "c", "map": {"52": 361, "53": 360, "56": 364, "64": 363}}]}]}
    partners = [{"id": "a", "ranges": [[348, 354, "H", "EEFFDAH"]]},
                {"id": "b", "ranges": [[343, 349, "H", "EEFFDAH"]]},
                {"id": "c", "ranges": [[358, 364, "H", "EEFFDAH"]]}]
    al = align_motif(entry, {"elements": [0], "carriers": ["a", "b", "c"]}, partners)
    cells = {c["id"]: c["cells"] for c in al["carriers"]}
    assert al["anchor"] == 52
    assert cells["a"] == {"-1": "F", "0": "F", "1": "D", "2": "A", "3": "H"}
    assert cells["c"] == {"-1": "F", "0": "F", "2": "A", "3": "H"}
    assert cells["b"] == {"-1": "F", "0": "F", "2": "A", "3": "H"}           # anchor alone put b at 0..+4


def test_logo_shows_every_touching_residue_not_only_the_mapped_one():
    """VAPA FFAT (2026-10-07): the engine keeps ONE partner residue per bait position (a: 348 on bait 54), so the
    motif's own E (349) vanished from the logo although it touches bait 54 too. `touching` adds it to the columns;
    placement is unchanged (it uses the map)."""
    entry = {"elements": [{"orientation": "parallel", "carriers": [
        {"id": "a", "group": "a", "map": {"52": 351, "54": 348}},
        {"id": "b", "group": "b", "map": {"52": 361, "54": 358}}]}]}
    partners = [{"id": "a", "ranges": [[348, 352, "H", "EEFFD"]]}, {"id": "b", "ranges": [[358, 362, "H", "EEFFD"]]}]
    motif = {"elements": [0], "carriers": ["a", "b"]}
    before = {c["id"]: c["cells"] for c in align_motif(entry, motif, partners)["carriers"]}
    assert before["a"] == {"-3": "E", "0": "F"}
    al = align_motif(entry, dict(motif, touching={"a": [348, 349, 350, 380], "b": [358, 359]}), partners)
    # 380 touches too but lies outside a's mapped span (348..351): left out, so the logo does not widen
    cells = {c["id"]: c["cells"] for c in al["carriers"]}
    assert cells["a"] == {"-3": "E", "-2": "E", "-1": "F", "0": "F"} and cells["b"] == {"-3": "E", "-2": "E", "0": "F"}
    col = next(c for c in al["columns"] if c.get("offset") == -2)
    assert col["k"] == 2 and col["counts"] == [["E", 2]]


def test_binding_subgroups_split_by_geometry_only():
    """2026-10-08 (viz_int 97_): carriers that touch the same bait residues with the same aligned residues form a
    binding sub-group. Six carriers touch bait 10/12 one way (a1-a3) or bait 14/16 another (b1-b3): two sub-groups.
    The letters are identical for all six, so only geometry can have split them."""
    from scrims.motif_align import binding_subgroups
    A = {"10": 100, "12": 102}; B = {"14": 100, "16": 102}
    cars = [{"id": f"{k}{i}", "group": f"{k}{i}", "map": dict(mp if k == "a" else B)}
            for k, mp in (("a", A), ("b", B)) for i in (1, 2, 3)]
    for c in cars:                      # everyone also touches bait 12 so the motif aligns as one
        c["map"].setdefault("12", 102)
    entry = {"elements": [{"orientation": "parallel", "carriers": cars}]}
    touching = {c["id"]: ({"10": [100], "12": [102]} if c["id"][0] == "a" else {"14": [100], "16": [102], "12": [102]})
                for c in cars}
    motif = {"elements": [0], "carriers": [c["id"] for c in cars], "touching": touching}
    al = align_motif(entry, motif, None)
    assert sorted(sorted(g) for g in al["subgroups"]) == [["a1", "a2", "a3"], ["b1", "b2", "b3"]]
    same = dict(motif, touching={c["id"]: {"10": [100], "12": [102]} for c in cars})
    assert binding_subgroups(align_motif(entry, same, None), same) == []    # one way of binding: no split


def test_motif_quality_counts_related_partners_once():
    """2026-10-09 (user option C): quality = bits of the 5 most conserved columns, paralog groups counted once.
    Three groups all with L at both columns: each column = 4.32 - e(3) bits. Adding two more members of group a
    (same letters) must not raise it -- a family repeated is one piece of evidence."""
    import math
    from scrims.logo import small_sample_correction
    from scrims.motif_align import motif_quality
    def run(extra):
        cars = [{"id": x, "group": x, "map": {"10": 100, "12": 102}} for x in ("a", "b", "c")]
        cars += [{"id": f"a{k}", "group": "a", "map": {"10": 100, "12": 102}} for k in range(extra)]
        parts = [{"id": c["id"], "ranges": [[100, 102, "H", "LxL"]]} for c in cars]
        al = align_motif({"elements": [{"orientation": "parallel", "carriers": cars}]},
                         {"elements": [0], "carriers": [c["id"] for c in cars]}, parts)
        return motif_quality(al)
    want = 2 * (math.log2(20) - small_sample_correction(3))
    assert abs(run(0) - want) < 1e-3 and abs(run(2) - want) < 1e-3

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


def test_quality_rank_against_the_fixed_reference():
    """'top N %' counts the reference motifs at least as informative, this one included (2026-10-12)."""
    from scrims.motif_align import QUALITY_REFERENCE, quality_top
    assert len(QUALITY_REFERENCE) == 51 and QUALITY_REFERENCE == sorted(QUALITY_REFERENCE)
    assert quality_top(99) == 2 and quality_top(0) == 100
    assert quality_top(6.699) == 27          # VAPA I1 FFAT: 13 reference motifs above it
