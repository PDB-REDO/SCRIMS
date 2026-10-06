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
  4. COLUMNS: every distinct partner residue a carrier has in the motif lands at
     dir x (residue - its reference residue) + reference offset. k = carriers with
     a residue there, n = aligned carriers. Runs of > MAX_EMPTY empty offsets
     collapse into one gap column.
Letters come from the partner annotation's `ranges`; a residue outside them is '?'.

Iteration order matters for step 3's tie-break and mirrors JavaScript: a map's
bait positions are visited in ascending order (JS orders integer-like object
keys that way), elements in the motif's order, carriers in first-seen order.
"""

from __future__ import annotations

MAX_EMPTY = 4


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

    lookup = _letter_lookup(partners)
    missing = 0
    by_offset, k_at, carriers = {}, {}, []
    for i in order:
        m = maps[i]
        residues = sorted(set(m.values()))
        seq = lookup.get(i)
        letters = ""
        for r in residues:
            ch = seq.get(r) if seq else None
            if not ch:
                missing += 1
            letters += ch or "?"
        r0, d = ref.get(i), dirs[i]
        if r0:
            for r, ch in zip(residues, letters):
                o = d * (r - r0[0]) + r0[1]
                cm = by_offset.setdefault(o, {})
                cm[ch] = cm.get(ch, 0) + 1
                k_at[o] = k_at.get(o, 0) + 1
        carriers.append({"id": i, **meta[i], "residues": residues, "letters": letters,
                         "aligned": r0 is not None, "dir": d})

    offs = sorted(by_offset)
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

    return {"anchor": anchor, "anchor_carriers": best,
            "orientation": "parallel" if fixed == 1 else "antiparallel" if fixed == -1 else "per-carrier",
            "n": sum(1 for c in carriers if c["aligned"]), "carriers": carriers, "columns": columns,
            "missing_letters": missing}
