"""The known-motif suggestion (scrims.motif_candidates): an ELM class the partners match beyond chance."""
import random

from scrims.motif_candidates import domain_evidence, motif_candidates, parse_domains
from scrims.structures import T2O

O2T = {v: k for k, v in T2O.items()}
AA = "ACDEFGHIKLMNPQRSTVWY"


def _screen(core, n=8, start=50, seed=0):
    """n unrelated partners, each a random 120-residue chain with `core` at `start`, mapped onto baits 10.."""
    rng = random.Random(seed)
    int_seq, maps = {}, {}
    for i in range(n):
        s = [rng.choice(AA) for _ in range(120)]
        if core:
            s[start:start + len(core)] = core
        int_seq[f"p{i}"] = {r + 1: O2T[a] for r, a in enumerate(s)}
        maps[f"p{i}"] = {10 + k: start + 1 + k for k in range(len(core) or 8)}
    elements = [{"map": maps}]
    motif = {"id": "M1", "positions": list(range(10, 18)), "elements": [0]}
    return motif, elements, int_seq, {c: c for c in int_seq}


def test_ffat_partners_get_the_ffat_class_as_likely():
    motif, elements, int_seq, groups = _screen("EEFFDAKE")
    got = motif_candidates("VAPA", "I1", motif, elements, {}, int_seq, groups)
    assert got and got[0]["elm"] == "TRG_ER_FFAT_1"
    c = got[0]
    assert c["groups_matched"] == 8 and c["groups"] == 8 and c["tier"] == "likely" and c["p"] <= 1e-6
    assert set(c["matches"]) == set(int_seq) and c["domain"] == "unknown"


def test_random_partners_get_no_suggestion():
    motif, elements, int_seq, groups = _screen("")
    assert motif_candidates("X", "I1", motif, elements, {}, int_seq, groups) == []


def test_related_partners_count_once():
    """Eight copies in ONE partner group are one independent match: no class can reach 3 groups."""
    motif, elements, int_seq, _ = _screen("EEFFDAKE")
    assert motif_candidates("VAPA", "I1", motif, elements, {}, int_seq, {c: "g" for c in int_seq}) == []


def test_domain_is_evidence_never_a_filter():
    motif, elements, int_seq, groups = _screen("EEFFDAKE")
    for doms, want in ((set(), "unknown"), (parse_domains("PF00635"), "supported"),
                       (parse_domains("IPR000535"), "unknown"), (parse_domains("PF00069"), "not found")):
        got = motif_candidates("VAPA", "I1", motif, elements, {}, int_seq, groups, doms)
        assert got[0]["elm"] == "TRG_ER_FFAT_1" and got[0]["domain"] == want, (doms, got[0]["domain"])


def test_parse_domains():
    assert parse_domains("PF00635.25, IPR000535 PF00069") == {"PF00635", "IPR000535", "PF00069"}
    assert parse_domains("") == set() and domain_evidence(set(), "TRG_ER_FFAT_1") == "unknown"
