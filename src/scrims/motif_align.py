"""
A motif's logo ALONG THE PARTNER SEQUENCE: the carriers of one motif aligned on
their own residues. Display only -- the engine never aligns partners; this is
the companion website's reading of the engine's carrier maps (its
`src/lib/motifAlign.ts`, 2026-10-02), ported line for line so the figures, the
HTML report and the website show the same logo.

Input per carrier: `map` {bait position -> partner residue}, from every element
of the motif it carries (merged; the first element in the motif's order wins a
clash on one bait position).
  1. ANCHOR = the motif's bait position held by the most carriers (ties -> the
     lowest position).
  2. DIRECTION: if all the motif's elements share one orientation (parallel /
     antiparallel) it holds for every carrier (antiparallel -> the partner is
     read C->N so it runs along the bait). Otherwise each carrier's own map
     decides: the sign of the sum over bait-position pairs of sign(partner-residue
     difference); a tie reads N->C.
  3. REGISTER: a carrier holding the anchor puts its anchor residue at offset 0.
     One that does not is placed through another bait position p it holds: at the
     lower-median offset that p's residue has among the anchor holders that also
     hold p (the p with the most such holders, first in iteration order on ties).
     A carrier with no such p is listed but left out of the logo (`aligned`
     False) -- never guessed.
  3b. CONSENSUS (2026-10-07, website and here alike): step 3 places each row through
     ONE bait position, which shifts a row whenever two equivalent partner residues
     swap on it (VAPA FFAT: both F of E-F-F-D-A touch bait 52 and 53; PITPNM2 was one
     column off PITPNM1 though offset -5 at 10 of its 11 positions). So: every bait
     position p held by an anchor holder gets a consensus offset c_p = the most common
     offset of p's residue among the holders (ties -> smallest |offset|, then the
     lower). Each carrier is then re-placed by a vote: every p it holds with a c_p
     proposes (direction, offset of its reference residue) that would put p's residue
     on c_p; the proposal with the most votes wins (ties -> keep the direction, then
     keep the offset, then the smallest change, then the lower offset). The direction
     is voted for every carrier (2026-10-07: also when all elements agree -- VAPA I3
     OSBPL7 binds antiparallel inside a "parallel" element and was forced forwards).
     Geometry only -- letters never place a row, so a logo cannot look conserved by
     construction. Measured (viz_int validation/scripts/95_register_consensus.py).
  3c. SIGNATURE PLACEMENT (2026-10-09, viz_int validation 98_): steps 3-3b read the engine's
     ONE-residue-per-bait-position map, which misplaces carriers whose map recorded another of
     several residues touching the anchor (PCNA: Q, F7 and more touch bait 252; POLD3 and RNASEH2B
     sat one column off, UNG reversed, although every canonical PIP-box touches PCNA identically).
     With motif["touching"] keyed by bait position, each carrier's (direction, offset) is re-chosen
     to maximise how many of its touching (bait, residue) pairs land where the OTHER carriers put
     theirs (consensus counts, own pairs left out); 3 rounds from the 3b placement, each round from
     the previous one; ties -> keep the direction, then the offset, then the smallest change, then
     the lower offset. Geometry only. Known motifs (20 classes): register 0.57 -> 0.65.
  3d. N->C DISPLAY (2026-10-09, user: "always show the N to C orientation ... keep the flagging"):
     if most independent groups read backwards, the whole motif is mirrored so it reads N->C;
     `flipped` = True flags that its partners run against the bait (antiparallel). Carriers read
     against the motif keep dir = -1 (shown <-). `subgroup_flip[k]` says whether binding sub-group k
     needs the same mirror when shown alone (its own majority).
  4. COLUMNS: every distinct partner residue a carrier has in the motif -- its map residues plus
     (2026-10-07) those of motif["touching"] (EVERY residue touching the motif's bait positions,
     since the engine keeps one per position and a conserved neighbour would vanish -- VAPA FFAT
     E348/E349 on bait 54) that lie between the carrier's first and last map residue -- lands at
     dir x (residue - its reference residue) + reference offset. k = carriers with
     a residue there, n = aligned carriers. Runs of > MAX_EMPTY empty offsets
     collapse into one gap column. Only offsets held by >= MIN_COL_GROUPS independent
     groups become columns (2026-10-07: one group = 0 bits, i.e. blank space).
Letters come from the partner annotation's `ranges`; a residue outside them is '?'.

Per carrier the result also keeps (website 2026-10-07) `map` {bait position: partner
residue} (merged over the motif's elements; drives the bait-side view filtered to one
motif) and `cells` {logo offset: the carrier's letter there} (drives the per-carrier
table and the group-weighted classic logo; empty for an unaligned carrier). Keys are
strings, as in JSON.

Iteration order matters for step 3's tie-break and mirrors JavaScript: a map's
bait positions are visited in ascending order (JS orders integer-like object
keys that way), elements in the motif's order, carriers in first-seen order.
"""

from __future__ import annotations

MAX_EMPTY = 4
MIN_COL_GROUPS = 2
SUBGROUP_CUT = 0.7


def _touching_of(motif, i):
    """A carrier's touching residues, flat (keyed by bait position since 2026-10-08; older documents: a list)."""
    t = (motif.get("touching") or {}).get(i)
    if not t:
        return []
    return list(t) if isinstance(t, list) else sorted({r for rs in t.values() for r in rs})


def _sign(x):
    return (x > 0) - (x < 0)


def _letter_lookup(partners):
    out = {}
    for p in partners or []:
        m = {}
        for a, b, _tier, seq in p["ranges"]:
            for r in range(a, b + 1):
                m[r] = seq[r - a] if r - a < len(seq) else "?"
        out[p.get("id", p.get("prediction"))] = m
    return out


def align_motif(entry, motif, partners):
    """`entry` = one interface's motif record ({'elements': [{'orientation', 'carriers':
    [{'id', 'group', 'map': {bait position: partner residue}}]}], ...}); `motif` = one of
    its motifs ({'elements': [element indices], 'carriers': [ids]}); `partners` = that
    interface's partner annotation (dicts with 'id' or 'prediction', and 'ranges').
    Returns the alignment dict: anchor, anchor_carriers, orientation, n, carriers, columns,
    missing_letters."""
    els = [entry["elements"][i] for i in motif["elements"] if 0 <= i < len(entry["elements"])]
    meta, maps = {}, {}
    for e in els:
        for c in e["carriers"]:
            meta.setdefault(c["id"], {"group": c.get("group"), "gene": c.get("gene"),
                                      "uniprot": c.get("uniprot")})
            m = maps.setdefault(c["id"], {})
            for p in sorted(c["map"], key=lambda k: int(k)):
                if int(p) not in m:
                    m[int(p)] = int(c["map"][p])
    order = [i for i in motif["carriers"] if i in maps] + [i for i in maps if i not in motif["carriers"]]

    ors = {e["orientation"] for e in els}
    fixed = 1 if ors == {"parallel"} else -1 if ors == {"antiparallel"} else None

    def dir_of(m):
        if fixed:
            return fixed
        ks = sorted(m)
        s = sum(_sign(m[ks[j]] - m[ks[i]]) for i in range(len(ks)) for j in range(i + 1, len(ks)))
        return 1 if s >= 0 else -1

    held = {}
    for m in maps.values():
        for p in m:
            held[p] = held.get(p, 0) + 1
    anchor, best = float("inf"), -1
    for p, n in held.items():
        if n > best or (n == best and p < anchor):
            anchor, best = p, n

    dirs = {i: dir_of(m) for i, m in maps.items()}
    holders = [i for i in order if anchor in maps[i]]
    ref = {i: (maps[i][anchor], 0) for i in holders}
    for i in order:
        if i in ref:
            continue
        pick = None                     # (n refs, partner residue, offset)
        for p, r in maps[i].items():
            offs = sorted(dirs[h] * (maps[h][p] - maps[h][anchor]) for h in holders if p in maps[h])
            if offs and (pick is None or len(offs) > pick[0]):
                pick = (len(offs), r, offs[(len(offs) - 1) // 2])
        if pick:
            ref[i] = (pick[1], pick[2])

    col = {}                            # 3b: consensus offset per bait position, from the anchor holders
    for h in holders:
        r0, o0 = ref[h]
        for p, r in maps[h].items():
            col.setdefault(p, {})
            o = dirs[h] * (r - r0) + o0
            col[p][o] = col[p].get(o, 0) + 1
    cons = {p: min(c, key=lambda o: (-c[o], abs(o), o)) for p, c in col.items()}
    for i in order:
        m = maps[i]
        r0, o0 = ref.get(i, (min(m.values()), None))
        votes = {}
        for p, r in m.items():
            if p not in cons:
                continue
            for d in (dirs[i], -dirs[i]):
                k = (d, cons[p] - d * (r - r0))
                votes[k] = votes.get(k, 0) + 1
        if votes:
            d, o = max(votes, key=lambda k: (votes[k], k[0] == dirs[i], k[1] == o0,
                                            -abs(k[1] - (o0 or 0)), -k[1]))
            dirs[i], ref[i] = d, (r0, o)

    # 3c: placement by touching signature (needs `touching` keyed by bait position)
    tp = {}
    for i in order:
        t = (motif.get("touching") or {}).get(i)
        if isinstance(t, dict) and t and i in ref:
            tp[i] = sorted({(int(b), int(r)) for b, rs in t.items() for r in rs})
    if len(tp) >= 2:
        place = {i: (dirs[i], ref[i][0], ref[i][1]) for i in tp}
        for _ in range(3):
            sig = {i: [(b, place[i][0] * (r - place[i][1]) + place[i][2]) for b, r in tp[i]] for i in tp}
            C = {}
            for i in tp:
                for bo in sig[i]:
                    C[bo] = C.get(bo, 0) + 1
            new = {}
            for i in tp:
                own = {}
                for bo in sig[i]:
                    own[bo] = own.get(bo, 0) + 1
                byb = {}
                for (b, o), n in C.items():
                    n -= own.get((b, o), 0)
                    if n > 0:
                        byb.setdefault(b, []).append((o, n))
                d0, r0, o0 = place[i]
                best = None
                for d in (d0, -d0):
                    votes = {}
                    for b, r in tp[i]:
                        for o, n in byb.get(b, ()):
                            sft = o - d * (r - r0)
                            votes[sft] = votes.get(sft, 0) + n
                    for sft, v in votes.items():
                        k = (v, d == d0, sft == o0, -abs(sft - o0), -sft)
                        if best is None or k > best[0]:
                            best = (k, (d, r0, sft))
                new[i] = best[1] if best else place[i]
            if new == place:
                break
            place = new
        for i, (d, r0, o0) in place.items():
            dirs[i], ref[i] = d, (r0, o0)

    # 3d: show the motif N->C; flag when its partners run against the bait
    def backwards(ids):
        net = {}
        for i in ids:
            net[meta[i]["group"]] = net.get(meta[i]["group"], 0) + dirs[i]
        return sum(v < 0 for v in net.values()) > sum(v > 0 for v in net.values())
    flipped = backwards([i for i in order if i in ref])
    if flipped:
        for i in ref:
            dirs[i], ref[i] = -dirs[i], (ref[i][0], -ref[i][1])

    lookup = _letter_lookup(partners)
    missing = 0
    by_offset, k_at, g_at, carriers = {}, {}, {}, []
    for i in order:
        m = maps[i]
        # placement (steps 1-3b) uses the map only; the columns add every touching residue (motif["touching"]) BETWEEN
        # the carrier's first and last mapped residue (fills the motif's gaps without widening the logo)
        lo, hi = min(m.values()), max(m.values())
        residues = sorted(set(m.values()) | {r for r in _touching_of(motif, i) if lo <= r <= hi})
        seq = lookup.get(i)
        letters = ""
        for r in residues:
            ch = seq.get(r) if seq else None
            if not ch:
                missing += 1
            letters += ch or "?"
        r0, d = ref.get(i), dirs[i]
        cells = {}
        if r0:
            for r, ch in zip(residues, letters):
                o = d * (r - r0[0]) + r0[1]
                cells[str(o)] = ch
                cm = by_offset.setdefault(o, {})
                cm[ch] = cm.get(ch, 0) + 1
                k_at[o] = k_at.get(o, 0) + 1
                g_at.setdefault(o, set()).add(meta[i]["group"])
        carriers.append({"id": i, **meta[i], "residues": residues, "letters": letters,
                         "aligned": r0 is not None, "dir": d,
                         "map": {str(p): r for p, r in sorted(m.items())}, "cells": cells,
                         "ref": list(r0) if r0 else None})

    offs = sorted(o for o in by_offset if len(g_at[o]) >= MIN_COL_GROUPS) or sorted(by_offset)
    columns = []
    for n, o in enumerate(offs):
        if n > 0:
            gap = o - offs[n - 1] - 1
            if gap > MAX_EMPTY:
                columns.append({"kind": "gap", "skipped": gap})
            else:
                columns += [{"kind": "col", "offset": offs[n - 1] + g, "k": 0, "counts": []}
                            for g in range(1, gap + 1)]
        counts = sorted(by_offset[o].items(), key=lambda kv: -kv[1])   # stable, as JS sort
        columns.append({"kind": "col", "offset": o, "k": k_at[o], "counts": [list(c) for c in counts]})

    out = {"anchor": anchor, "anchor_carriers": best,
            "orientation": "parallel" if fixed == 1 else "antiparallel" if fixed == -1 else "per-carrier",
            "n": sum(1 for c in carriers if c["aligned"]), "carriers": carriers, "columns": columns,
            "missing_letters": missing, "flipped": flipped}
    out["subgroups"] = binding_subgroups(out, motif)
    out["subgroup_flip"] = [backwards(g) for g in out["subgroups"]]
    out["quality"] = motif_quality(out)
    out["quality_top"] = quality_top(out["quality"])
    out["consensus"] = motif_consensus(out)
    out["best_example"] = best_example(out)
    # per binding sub-group, in the orientation its view is shown in (mirrored when subgroup_flip)
    out["subgroup_consensus"], out["subgroup_best_example"] = [], []
    for g, fl in zip(out["subgroups"], out["subgroup_flip"]):
        v = mirror(out) if fl else out
        out["subgroup_consensus"].append(motif_consensus(v, set(g)))
        out["subgroup_best_example"].append(best_example(v, set(g)))
    return out


def _col_freq(cars, offset):
    from scrims.logo import group_weighted
    return group_weighted([(c["cells"][str(offset)], c["group"]) for c in cars if str(offset) in c["cells"]])


SIMILAR = set("AS DE DN EK EQ FW FY HN HY IL IM IV KQ KR LM LV MV NS QR ST WY".split())   # amino-acid pairs with BLOSUM62 >= 1 (Biopython BLOSUM62)


def _bracket(top):
    import itertools
    pool = [(l, v) for l, v in top[:6] if v >= 0.1]
    for k in (2, 3):
        best = None
        for combo in itertools.combinations(range(len(pool)), k):
            ls = [pool[i][0] for i in combo]
            if not all((a + b if a < b else b + a) in SIMILAR for a, b in itertools.combinations(ls, 2)):
                continue
            s = sum(pool[i][1] for i in combo)
            if s >= 0.6 - 1e-12 and (best is None or s > best[1] + 1e-12):
                best = (ls, s)
        if best:
            return "[" + "".join(best[0]) + "]"
    return None


def motif_consensus(al, ids=None):
    """CONSENSUS (2026-10-09, user option A; website motifAlign.ts motifConsensus, identical): [[offset | None, token]],
    counted like the classic logo (every carrier). Core columns = held by >= max(2, 25 % of the carriers); per core
    column: the top letter when >= 50 %; else the smallest set of 2 (then 3) letters, each >= 10 %, together >= 60 %,
    pairwise BLOSUM62 >= 1 (the engine's similarity rule), in brackets; else x. Non-core columns inside the core span
    are x, a collapsed gap is …; leading/trailing x trimmed."""
    import math
    cars = [c for c in al["carriers"] if c["aligned"] and (ids is None or c["id"] in ids)]
    need = max(2, math.ceil(0.25 * len(cars)))
    toks = []
    for col in al["columns"]:
        if col["kind"] == "gap":
            toks.append([None, "…"]); continue
        cnt = {}
        for c in cars:
            l = c["cells"].get(str(col["offset"]))
            if l:
                cnt[l] = cnt.get(l, 0) + 1
        n = sum(cnt.values())
        top = sorted(((l, k / max(1, n)) for l, k in cnt.items()), key=lambda t: (-t[1], t[0]))
        tok = ""
        if n >= need:
            tok = top[0][0] if top[0][1] >= 0.5 else (_bracket(top) or "x")
        toks.append([col["offset"], tok])
    core = [i for i, t in enumerate(toks) if t[1] not in ("", "x", "…")]
    if not core:
        return []
    return [[o, t or "x"] for o, t in toks[core[0]:core[-1] + 1]]


def best_example(al, ids=None):
    """BEST EXAMPLE (2026-10-09, user option C; website bestExample, identical): the carrier id whose letters agree most
    with the consensus columns (sum of the group-weighted share of its letter); ties -> more columns covered, then id."""
    cars = [c for c in al["carriers"] if c["aligned"] and (ids is None or c["id"] in ids)]
    cons = [o for o, t in motif_consensus(al, ids) if o is not None and t not in ("x", "…")]
    if not cons:
        return None
    share = {o: dict(_col_freq(cars, o)[0]) for o in cons}
    best = None
    for c in cars:
        sc = cov = 0
        for o in cons:
            l = c["cells"].get(str(o))
            if l:
                sc += share[o].get(l, 0.0); cov += 1
        k = (round(sc, 9), cov, c["id"])
        if best is None or k[0] > best[0] or (k[0] == best[0] and (k[1] > best[1] or (k[1] == best[1] and k[2] < best[2]))):
            best = k
    return best[2]


QUALITY_TOP = 5


def motif_quality(al, ids=None):
    """MOTIF QUALITY (2026-10-09, user option C): the information (bits) of the motif's QUALITY_TOP most conserved logo
    columns, related partners (paralog groups) counted once -- so neither a long stretch (many weak columns) nor a
    protein family (one sequence repeated) inflates it. Per column: group-weighted frequencies and the exact
    small-sample correction (scrims.logo, the bait-side classic logo's maths); only columns held by >= 2 groups.
    Over the 31-protein bundle: 0.6-10.7, median 4.3. Same as website motifAlign.ts motifQuality."""
    from scrims.logo import group_weighted, info_content
    cars = [c for c in al["carriers"] if c["aligned"] and (ids is None or c["id"] in ids)]
    R = []
    for col in al["columns"]:
        if col["kind"] != "col":
            continue
        items = [(c["cells"][str(col["offset"])], c["group"]) for c in cars if str(col["offset"]) in c["cells"]]
        freq, n = group_weighted(items)
        if n >= 2:
            R.append(info_content([f for _, f in freq], n))
    return round(sum(sorted(R, reverse=True)[:QUALITY_TOP]), 3)


# QUALITY RANK (2026-10-12, user's PI: "add significance to the bits"; user: options A + C): the bits are a size, not
# a test, so the report shows them with the motif's rank among a FIXED reference set -- the quality of the 51 motifs of
# viz_int's 31-protein bundle (website build 5QTT16Cxi5JSHNQT6XX_r, 2026-10-12) -- and the significance the motif already
# has (its strongest element's BH q, report/motifs.tsv). Same as website motifAlign.ts QUALITY_REFERENCE / qualityTop.
QUALITY_REFERENCE = [0.568, 0.568, 1.135, 1.395, 1.499, 1.518, 1.895, 1.895, 2.054, 2.054, 2.084, 2.253, 2.271, 2.272, 2.377, 2.621, 2.621, 2.621, 2.651, 2.972, 3.189, 3.189, 3.54, 3.844, 4.128, 4.295, 4.373, 4.514, 4.532, 4.656, 5.198, 5.419, 5.44, 6.206, 6.304, 6.505, 6.596, 6.699, 6.719, 6.894, 7.261, 7.602, 7.72, 7.957, 8.245, 8.433, 8.867, 8.889, 9.279, 10.08, 10.731]


def quality_top(q):
    """'top N %': the share of reference motifs at least as informative, counting this one (1..100)."""
    above = sum(1 for x in QUALITY_REFERENCE if x > q + 1e-9)
    return max(1, round(100 * (above + 1) / (len(QUALITY_REFERENCE) + 1)))


def mirror(al):
    """The same alignment read the other way round (offsets negated, directions reversed): what a binding sub-group
    whose own majority reads backwards is shown with (`subgroup_flip`)."""
    cars = [dict(c, dir=-c["dir"], cells={str(-int(o)): l for o, l in c["cells"].items()},
                 ref=[c["ref"][0], -c["ref"][1]] if c.get("ref") else None) for c in al["carriers"]]
    cols = [dict(c, offset=-c["offset"]) if c["kind"] == "col" else c for c in reversed(al["columns"])]
    return dict(al, carriers=cars, columns=cols, flipped=not al.get("flipped", False))


def binding_subgroups(al, motif):
    """BINDING SUB-GROUPS (viz_int 2026-10-08, validation 97_; website motifAlign.ts bindingSubgroups, identical):
    carriers that touch the same bait residues with the same aligned residues bind the same way. Signature = every
    (bait position, logo offset) of a carrier's touching residues (its map when `touching` is absent or flat);
    distance = 1 - Jaccard; average linkage, merged while the closest pair of clusters is <= SUBGROUP_CUT apart (ties ->
    lowest indices). Geometry only -- letters never decide. Returns [[carrier ids]] for clusters with >= 3 independent
    groups, largest first; [] when the only such cluster holds >= 75 % of the aligned carriers (no split)."""
    cars = [c for c in al["carriers"] if c["aligned"] and c.get("ref")]
    sig = []
    for c in cars:
        r0, o0 = c["ref"]; t = (motif.get("touching") or {}).get(c["id"])
        pairs = ([(b, r) for b, rs in t.items() for r in rs] if isinstance(t, dict) and t
                 else [(b, r) for b, r in c["map"].items()])
        sig.append({f"{int(b)}:{c['dir'] * (r - r0) + o0}" for b, r in pairs})
    n = len(cars)
    D = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            inter = len(sig[i] & sig[j]); uni = len(sig[i]) + len(sig[j]) - inter
            D[i][j] = D[j][i] = 1 - (inter / uni if uni else 0)
    clusters = [[i] for i in range(n)]
    while True:
        bi = bj = -1; bd = float("inf")
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                d = sum(D[x][y] for x in clusters[i] for y in clusters[j]) / (len(clusters[i]) * len(clusters[j]))
                if d < bd - 1e-12:
                    bd, bi, bj = d, i, j
        if bi < 0 or bd > SUBGROUP_CUT + 1e-12:
            break
        merged = sorted(clusters[bi] + clusters[bj])
        clusters = [c for k, c in enumerate(clusters) if k not in (bi, bj)] + [merged]
    big = sorted((c for c in clusters if len({cars[i]["group"] for i in c}) >= 3), key=lambda c: (-len(c), c[0]))
    if len(big) == 1 and len(big[0]) >= 0.75 * n:
        return []
    return [[cars[i]["id"] for i in c] for c in big]
