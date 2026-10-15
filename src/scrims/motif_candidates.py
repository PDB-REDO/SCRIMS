"""
Known-class suggestion for a discovered motif: which ELM class (if any) its partner sequences match beyond chance.

Port of viz_int.motif_candidates (2026-10-10; viz_int IMPROVEMENTS #41, validated in viz_int
validation/scripts/99_motif_candidates.py and 99b_candidates_shipped.py) with one change: the bait's domain
annotation comes from the run's `bait_domains` parameter (Pfam PF... / InterPro IPR... accessions) instead of
viz_int's per-gene tables -- SCRIMS has no annotation of its own, so the evidence is "unknown" unless given.

A motif is suggested an ELM class when its carriers match that class's pattern beyond chance:
  - carrier stretch = its motif residues (element maps, plus touching residues between its first and last mapped one --
    the residues the motif panel shows) +- FLANK, from the full partner chain;
  - every ELM binding class (data/elms_index.tsv, ELM 1.4 download 2026-09-29; MOD_ and CLV_ excluded);
  - a GROUP (related partners) matches if any member's stretch matches; chance per member = (k + 1) / (N_SHUF + 1), k =
    composition-preserving shuffles of its own stretch that match (seeded per bait/interface/motif/carrier; the +1
    keeps an estimate from 200 shuffles from claiming p = 0); per group 1 - prod(1 - p_member); observed matched groups
    vs the Poisson-binomial of those chances -> P(X >= observed);
  - suggested when p < 0.05 / (classes tested in this motif) (Bonferroni within the motif), observed >= 2 x expected,
    and >= 3 groups or >= half of them match;
  - tier "likely" when >= 5 groups match and p <= 1e-6, else "possible" (few carriers cannot separate a rare real motif
    from a lucky match -- the tier says so);
  - domain evidence (never a filter): "supported" when the bait carries a domain ELM lists for the class (same Pfam
    accession or clan, or InterPro id), "not found" when the bait is annotated in the scheme ELM uses and lacks it,
    "unknown" otherwise.
Measured in viz_int (31 proteins): known class named first in 9 of 10 motifs holding one; 12 of 41 other motifs get a
suggestion; the "likely" tier (5) all correct or plausible; every chance hit carries "not found".

Data: ELM classes and interaction domains (elm.eu.org, 1.4, 2026-09-29); Pfam-A clans (Pfam, CC0).
"""

from __future__ import annotations

import collections
import csv
import gzip
import hashlib
import random
import re
from functools import lru_cache
from pathlib import Path

import numpy as np

from scrims.structures import aa1

DATA = Path(__file__).parent / "data"
FLANK, N_SHUF, ALPHA, FOLD, TOP = 8, 200, 0.05, 2.0, 2
LIKELY_GROUPS, LIKELY_P = 5, 1e-6


@lru_cache(maxsize=1)
def _tables():
    rows = csv.DictReader((l for l in open(DATA / "elms_index.tsv") if not l.startswith("#")), delimiter="\t")
    elm = []
    for r in rows:
        if r["ELMIdentifier"].startswith(("MOD_", "CLV_")):
            continue
        try:
            elm.append((r["ELMIdentifier"], r["FunctionalSiteName"].strip(), re.compile(r["Regex"])))
        except re.error:
            pass
    dom = collections.defaultdict(set)
    for r in csv.DictReader((l for l in open(DATA / "elm_interaction_domains.tsv") if not l.startswith("#")), delimiter="\t"):
        dom[r["ELM identifier"]].add(r["Interaction Domain Id"].split(".")[0])
    clan = {}
    for line in gzip.open(DATA / "Pfam-A.clans.tsv.gz", "rt"):
        f = line.rstrip("\n").split("\t")
        clan[f[0]] = f[1] or f[0]
    return elm, dom, clan


def parse_domains(spec):
    """'PF00635,IPR000535' (comma/space separated; version suffixes dropped) -> {accession}."""
    return {a.split(".")[0] for a in re.split(r"[,\s]+", spec or "") if a}


def domain_evidence(bait_domains, elm_id):
    """'supported' / 'not found' / 'unknown': does the bait carry a domain ELM lists for this class?"""
    _, dom, clan = _tables()
    doms = dom.get(elm_id, set())
    pf = {a for a in doms if a.startswith("PF")}
    ipr = {a for a in doms if a.startswith("IPR")}
    bpf = {a for a in bait_domains if a.startswith("PF")}; bcl = {clan.get(a, a) for a in bpf}
    bip = {a for a in bait_domains if a.startswith("IPR")}
    if any(a in bpf or clan.get(a, a) in bcl for a in pf) or (ipr & bip):
        return "supported"
    if (pf and bpf) or (ipr and bip):
        return "not found"
    return "unknown"


def _pb_tail(ps, k):
    dist = np.zeros(len(ps) + 1); dist[0] = 1.0
    for p in ps:
        dist[1:] = dist[1:] * (1 - p) + dist[:-1] * p
        dist[0] *= (1 - p)
    return float(dist[k:].sum())


def motif_candidates(seed, iface, motif, elements, touch, int_seq, groups, bait_domains=frozenset()):
    """[{elm, name, regex, tier, groups_matched, groups, expected, p, domain, matches: {carrier: [start, end]}}] (lowest p
    first, <= TOP) for one engine motif. `seed` = the bait's name (gene, else the screen name; seeds the shuffles),
    `elements` = the interface's engine elements (map per carrier), `touch` = {carrier: {bait: [residues]}}
    (screen.py), `int_seq` = {carrier: {residue: three-letter name}}, `groups` = {carrier: group}."""
    elm, _, _ = _tables()
    pos = {int(x) for x in motif["positions"]}
    maps = collections.defaultdict(dict)
    for k in motif["elements"]:
        for c, m in elements[k]["map"].items():
            for b, r in m.items():
                maps[c].setdefault(int(b), int(r))
    seg, shuf = {}, {}
    for c, m in maps.items():
        seq = int_seq.get(c) or {}
        vals = list(m.values()); lo, hi = min(vals), max(vals)
        tr = {int(r) for b, rs in (touch.get(c) or {}).items() if int(b) in pos for r in rs if lo <= int(r) <= hi}
        res = set(vals) | tr
        rr = [r for r in range(min(res) - FLANK, max(res) + FLANK + 1) if r in seq]
        s = "".join(aa1(seq[r]) for r in rr)
        seg[c] = (s, rr)
        rng = random.Random(int(hashlib.sha1(f"{seed}|{iface}|{motif['id']}|{c}".encode()).hexdigest()[:12], 16))
        letters = list(s); shuf[c] = []
        for _ in range(N_SHUF):
            rng.shuffle(letters)
            shuf[c].append("".join(letters))
    byg = collections.defaultdict(list)
    for c in seg:
        byg[groups.get(c, c)].append(c)
    tested = []
    for eid, name, rx in elm:
        hit = {c: rx.search(seg[c][0]) for c in seg}
        obs = sum(any(hit[c] for c in mem) for mem in byg.values())
        if obs < 2:
            continue
        pg = []
        for mem in byg.values():
            q = 1.0
            for c in mem:
                q *= 1 - (sum(1 for s in shuf[c] if rx.search(s)) + 1) / (N_SHUF + 1)
            pg.append(1 - q)
        exp = sum(pg)
        matches = {c: [seg[c][1][h.start()], seg[c][1][h.end() - 1]] for c, h in hit.items() if h}
        tested.append(dict(elm=eid, name=name, regex=rx.pattern, groups_matched=obs, groups=len(byg), expected=round(exp, 2),
                           p=_pb_tail(pg, obs), matches=matches))
    cut = ALPHA / max(1, len(tested))
    good = [t for t in tested
            if t["p"] < cut and t["groups_matched"] >= FOLD * max(t["expected"], 1e-9)
            and (t["groups_matched"] >= 3 or t["groups_matched"] >= 0.5 * t["groups"])]
    good.sort(key=lambda t: (t["p"], -t["groups_matched"]))
    out = []
    for t in good[:TOP]:
        t["tier"] = "likely" if t["groups_matched"] >= LIKELY_GROUPS and t["p"] <= LIKELY_P else "possible"
        t["domain"] = domain_evidence(bait_domains, t["elm"])
        t["p"] = float(f"{t['p']:.3g}")
        out.append(t)
    return out
