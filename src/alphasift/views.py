"""
Static figures of the companion website's views, drawn from an `alphasift/1`
document (alphasift.export.document), so the figures, the JSON and the HTML report
show the same thing.

  overview(doc, path)            confidence matrix (binder rows, clustered; contacted
                                 residues as columns) + RMSF and recurrence tracks +
                                 one interface lane
  interface_card(doc, i, path)   one interface: bait-side view (Panel 1 letter columns,
                                 height ~ sqrt(groups), * q < 0.05, hatched = untested),
                                 motif lanes, the relations line, then per motif its
                                 prey-sequence logo, carriers (by paralog group, with
                                 reading direction) and element table
  partner_strips(doc, i, path)   every prey at the interface with its tiered residue
                                 ranges (H solid, M light, L dashed), motif residues bold

Colours, sizes and rules follow the website (webiste_viz_int src/lib/colors.ts and
src/components/{Heatmap,MotifEngineView}.tsx, 2026-10-02): Greens ramp for confidence,
purple ramp for RMSF, eight categorical interface hues keyed by interface order
(significant = solid, not = hatched), amino acids in Zappo groups with H separate.
Colour is never the only cue: every band, column and chip also carries text.

Layout is in pixels on one axes per figure (1 px = 1 data unit at dpi 100); `dpi`
sets the output resolution.
"""

from __future__ import annotations

import math

import numpy as np

from alphasift.export import label

CONF_LO, CONF_HI = "#f7fcf5", "#00441b"
RMSF_LO, RMSF_HI = "#efedf5", "#3f007d"
HIGHLIGHT = "#ffc400"
INTERFACE_COLORS = ["#2a78d6", "#d06f00", "#d4609a", "#a0662f", "#1a9fc0", "#d4462c", "#7b61d9", "#c2566b"]
AA_GROUPS = [("aliphatic / hydrophobic", "ILVAM", "#242c35"), ("aromatic", "FWY", "#a169eb"),
             ("positive", "KR", "#2d3ccb"), ("histidine", "H", "#0e9f93"), ("negative", "DE", "#ec5560"),
             ("hydrophilic", "STNQ", "#335309"), ("conformationally special", "PG", "#904680"),
             ("cysteine", "C", "#847f13")]
AA_COLOR = {a: c for _, letters, c in AA_GROUPS for a in letters}
INK, INK2, INK3, LINE = "#191c20", "#4a5058", "#7c828b", "#e2e4e8"
Q_SIG = 0.05
TIER_STYLE = {"H": dict(fc="#25603c", ec="#25603c", ls="-", tc="white"),     # the website's tier chips
              "M": dict(fc="#adc5b4", ec="#adc5b4", ls="-", tc=INK),
              "L": dict(fc="white", ec=INK3, ls="--", tc=INK2)}


def aa_color(ch):
    return AA_COLOR.get(ch.upper(), "#8a9099")


def iface_color(idx):
    return INTERFACE_COLORS[idx % len(INTERFACE_COLORS)]


def _lerp(lo, hi, t):
    a = np.array([int(lo[i:i + 2], 16) for i in (1, 3, 5)], float)
    b = np.array([int(hi[i:i + 2], 16) for i in (1, 3, 5)], float)
    return tuple((a + (b - a) * min(1.0, max(0.0, t))) / 255.0)


def fmt_q(q):
    return "< 0.001" if q < 0.001 else f"{q:.2g}"


def compact(rs):
    """'73, 75, 77-78' -- residue numbers, runs compacted."""
    rs, out, i = sorted(rs), [], 0
    while i < len(rs):
        j = i
        while j + 1 < len(rs) and rs[j + 1] == rs[j] + 1:
            j += 1
        out.append(f"{rs[i]}–{rs[j]}" if j > i else str(rs[i]))
        i = j + 1
    return ", ".join(out)


class _Canvas:
    """A figure whose axes is a W x H pixel grid with y pointing down."""

    def __init__(self, W, H):
        import matplotlib.pyplot as plt
        self.W, self.H = W, H
        self.fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, W)
        self.ax.set_ylim(H, 0)
        self.ax.axis("off")

    def rect(self, x, y, w, h, **kw):
        from matplotlib.patches import Rectangle
        self.ax.add_patch(Rectangle((x, y), w, h, **kw))

    def text(self, x, y, s, size=8, color=INK, ha="left", va="baseline", weight="normal", **kw):
        self.ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, fontweight=weight, **kw)

    def line(self, xs, ys, **kw):
        self.ax.plot(xs, ys, **kw)

    def glyph(self, ch, x, y_bottom, w, h, color):
        """Letter `ch` stretched into the w x h box whose bottom edge is at y_bottom
        (the website's Glyph: same horizontal scale for every letter, wide ones squeezed)."""
        if h < 1.2:
            return
        from matplotlib.font_manager import FontProperties
        from matplotlib.patches import PathPatch
        from matplotlib.textpath import TextPath
        from matplotlib.transforms import Affine2D
        fp = FontProperties(family="DejaVu Sans", weight="bold")
        p = TextPath((0, 0), ch, size=10, prop=fp)
        cap = TextPath((0, 0), "H", size=10, prop=fp).get_extents()
        bb = p.get_extents()
        sx = min(0.9 * w / max(cap.width, 1e-6), 0.95 * w / max(bb.width, 1e-6))
        sy = h / max(bb.height, 1e-6)
        t = (Affine2D().translate(-bb.x0 - bb.width / 2, -bb.y0).scale(sx, -sy)
             .translate(x + w / 2, y_bottom))
        self.ax.add_patch(PathPatch(t.transform_path(p), fc=color, ec="none", lw=0))

    def crop(self, h):
        self.H = int(h)
        self.fig.set_size_inches(self.W / 100, self.H / 100)
        self.ax.set_ylim(self.H, 0)

    def save(self, path, dpi=200):
        self.fig.savefig(path, dpi=dpi)
        import matplotlib.pyplot as plt
        plt.close(self.fig)


def _bait_name(doc):
    idn = doc["bait"].get("identity") or {}
    return idn.get("gene") or doc["screen"]["name"]


# ----------------------------------------------------------------------------- overview
def overview(doc, path, dpi=200):
    """The website's heatmap block: RMSF and recurrence tracks, the confidence matrix
    (binder rows in clustered order, one column per bait residue some binder contacts),
    and one interface lane under the axis."""
    sig = np.array(doc["matrix"]["signal"], float)
    rows = doc["matrix"]["row_order"]
    if not rows:
        raise ValueError("no binder rows: nothing to draw")
    M = sig[rows]
    contacted = np.where((M > 0).any(0))[0]
    residues = [doc["matrix"]["residues"][c] for c in contacted]
    H_ = M[:, contacted]
    n_b, n_c = H_.shape
    col_of = {r: i for i, r in enumerate(residues)}

    PAD_L, PAD_R, PAD_T, TRACK_H, GAP = 70, 16, 34, 13, 3
    W = 1000
    grid_w = W - PAD_L - PAD_R
    cw = grid_w / max(1, n_c)
    row_h = max(1.0, min(3.4, 340 / n_b))
    mat_top = PAD_T + 2 * TRACK_H + GAP + 8
    mat_h = n_b * row_h
    axis_y = mat_top + mat_h
    lane_y = axis_y + 40
    H = int(lane_y + 12 + 70)
    cv = _Canvas(W, H)

    name = _bait_name(doc)
    cv.text(PAD_L, 18, f"{name}  ·  {doc['bait']['length']} aa  ·  {n_b} binders of "
            f"{doc['screen']['n_predictions']} predictions  ·  {len(doc['interfaces'])} interfaces",
            size=10, weight="bold")

    # tracks: RMSF (normalised per protein) and recurrence = summed signal / its maximum
    rm = dict(zip(doc["tracks"]["rmsf"]["residues"], doc["tracks"]["rmsf"]["values"]))
    rmax = max(rm.values()) if rm else 1.0
    ssum = np.array(doc["tracks"]["summed_signal"], float)
    smax = max(1e-6, float(ssum[contacted].max())) if n_c else 1.0
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    for k, r in enumerate(residues):
        x = PAD_L + k * cw
        cv.rect(x, PAD_T, cw + 0.3, TRACK_H, fc=_lerp(RMSF_LO, RMSF_HI, rm.get(r, 0.0) / (rmax or 1)), lw=0)
        cv.rect(x, PAD_T + TRACK_H + GAP, cw + 0.3, TRACK_H,
                fc=_lerp(CONF_LO, CONF_HI, ssum[r - 1] / smax), lw=0)
    cv.text(PAD_L - 6, PAD_T + TRACK_H - 3, "RMSF", size=7, color=INK3, ha="right", family="monospace")
    cv.text(PAD_L - 6, PAD_T + 2 * TRACK_H + GAP - 3, "recur", size=7, color=INK3, ha="right", family="monospace")

    # matrix as an image (values <= 0.02 left blank, as on the website)
    img = np.ones((n_b, n_c, 3))
    lo = np.array(_lerp(CONF_LO, CONF_HI, 0)); hi = np.array(_lerp(CONF_LO, CONF_HI, 1))
    img = lo + (hi - lo) * np.clip(H_, 0, 1)[..., None]
    img[H_ <= 0.02] = 1.0
    cv.ax.imshow(img, extent=[PAD_L, PAD_L + grid_w, axis_y, mat_top], aspect="auto",
                 interpolation="nearest", zorder=1)
    cv.rect(PAD_L, mat_top, grid_w, mat_h, fc="none", ec=LINE, lw=0.6, zorder=2)
    cv.text(PAD_L - 30, mat_top + mat_h / 2, f"{n_b} binders", size=7, color=INK3, ha="center",
            va="center", rotation=90, family="monospace")
    step = max(1, round(n_c / 8))
    for k in range(0, n_c, step):
        cv.text(PAD_L + (k + 0.5) * cw, axis_y + 12, str(residues[k]), size=7, color=INK3, ha="center",
                family="monospace")
    cv.text(PAD_L + grid_w / 2, axis_y + 26, f"{name} residue — columns: residues some binder contacts",
            size=7, color=INK3, ha="center", family="monospace")

    # interface lane: every interface's segments on one row, in its colour; hatched = not significant
    cv.text(PAD_L - 6, lane_y + 9, "interfaces", size=7, color=INK3, ha="right", family="monospace")
    cv.rect(PAD_L, lane_y, grid_w, 12, fc="#f3f4f6", lw=0)
    for idx, i in enumerate(doc["interfaces"]):
        c = iface_color(idx)
        xs = []
        for a, b in i["segments"]:
            cols = [col_of[r] for r in range(a, b + 1) if r in col_of]
            if not cols:
                continue
            x0, x1 = PAD_L + min(cols) * cw, PAD_L + (max(cols) + 1) * cw
            xs.append(x0)
            cv.rect(x0, lane_y, max(2, x1 - x0), 12, fc=c if i["significant"] else c + "38",
                    ec=c, lw=0 if i["significant"] else 0.8, hatch=None if i["significant"] else "////")
        if xs:
            cv.text(min(xs), lane_y + 24, i["id"], size=7, color=c, weight="bold")

    # legend
    ly = H - 22
    x = PAD_L
    cv.text(x, ly + 8, "confidence", size=7, color=INK2)
    for k in range(40):
        cv.rect(x + 60 + k * 1.5, ly, 1.6, 9, fc=_lerp(CONF_LO, CONF_HI, k / 39), lw=0)
    cv.text(x + 124, ly + 8, f"low → high (0.5 = PAE {doc['params']['d0_signal']:g} Å)", size=7, color=INK2)
    x += 330
    cv.rect(x, ly, 14, 9, fc=INTERFACE_COLORS[0], lw=0)
    cv.text(x + 18, ly + 8, "significant", size=7, color=INK2)
    cv.rect(x + 80, ly, 14, 9, fc=INTERFACE_COLORS[0] + "38", ec=INTERFACE_COLORS[0], lw=0.8, hatch="////")
    cv.text(x + 98, ly + 8, "not significant", size=7, color=INK2)
    x += 200
    cv.text(x, ly + 8, "RMSF", size=7, color=INK2)
    for k in range(40):
        cv.rect(x + 34 + k * 1.5, ly, 1.6, 9, fc=_lerp(RMSF_LO, RMSF_HI, k / 39), lw=0)
    cv.text(x + 98, ly + 8, "low → high", size=7, color=INK2)
    cv.save(path, dpi)


# ----------------------------------------------------------------------- interface card
COL_W, JUMP_W, ROW_GROUPS, ROW_SIG, LOGO_H, ROW_BAIT, ROW_NUM, LANE_H, LANE_GAP = 16, 5, 12, 12, 80, 16, 26, 13, 4
PCOL_W, GAP_W, GUTTER = 18, 30, 110
MIN_SHARE = 0.04


def _bait_side_height(entry):
    return ROW_GROUPS + ROW_SIG + LOGO_H + ROW_BAIT + ROW_NUM + 4 + len(entry["motifs"]) * (LANE_H + LANE_GAP)


def _draw_bait_side(cv, entry, seq, x0, y0):
    cols = entry["panel1"]
    xs, x = [], 0
    for k, c in enumerate(cols):
        if k > 0 and c["pos"] - cols[k - 1]["pos"] > 1:
            x += JUMP_W
        xs.append(x)
        x += COL_W
    max_g = max([1] + [c["n_groups"] for c in cols])
    y_logo = y0 + ROW_GROUPS + ROW_SIG
    y_bait = y_logo + LOGO_H
    y_num = y_bait + ROW_BAIT
    y_lanes = y_num + ROW_NUM + 4
    for lab, yy, hh in (("groups", y0, ROW_GROUPS), (f"✱ q < {Q_SIG}", y0 + ROW_GROUPS, ROW_SIG),
                        ("partner residues", y_logo, LOGO_H), ("bait residue", y_bait, ROW_BAIT),
                        ("position", y_num, ROW_NUM)):
        cv.text(x0 - 8, yy + hh / 2, lab, size=6.5, color=INK3, ha="right", va="center")
    for c, cx in zip(cols, xs):
        X = x0 + cx
        untested = c["q"] is None
        sig = not untested and c["q"] < Q_SIG
        col_h = LOGO_H * math.sqrt(c["n_groups"] / max_g)
        if sig:
            cv.rect(X, y_logo - 1, COL_W, LOGO_H + 2, fc="#e6ebf2", lw=0)
            cv.rect(X + 1, y_logo + LOGO_H, COL_W - 2, 2, fc=INK, lw=0)
            cv.text(X + COL_W / 2, y0 + ROW_GROUPS + ROW_SIG - 2, "✱", size=7, color=INK, ha="center")
        if untested:
            cv.rect(X, y_logo, COL_W, LOGO_H, fc="none", ec="#d6d9de", lw=0, hatch="////")
        cv.text(X + COL_W / 2, y0 + ROW_GROUPS - 3, str(c["n_groups"]), size=6, ha="center",
                color="#a3a9b1" if untested else INK2)
        y = y_logo + LOGO_H
        for l, s in sorted(c["letters"].items(), key=lambda kv: kv[1]):     # smallest first, largest on top
            if s < MIN_SHARE:
                continue
            h = s * col_h
            cv.glyph(l, X, y, COL_W, h, "#b3b8bf" if untested else aa_color(l))
            y -= h
        b = seq[c["pos"] - 1] if seq and c["pos"] - 1 < len(seq) else "?"
        cv.text(X + COL_W / 2, y_bait + 12, b, size=8, color=aa_color(b), ha="center", weight="bold")
        cv.text(X + COL_W / 2, y_num + 2, str(c["pos"]), size=6, color=INK3, ha="center", va="top",
                rotation=90, family="monospace")
    col_of = {c["pos"]: k for k, c in enumerate(cols)}
    for li, m in enumerate(entry["motifs"]):
        y = y_lanes + li * (LANE_H + LANE_GAP)
        idx = [col_of[p] for p in m["positions"] if p in col_of]
        cv.text(x0 - 8, y + LANE_H / 2, m["id"], size=7, color=INK, ha="right", va="center", weight="bold")
        if not idx:
            continue
        a, b = x0 + xs[min(idx)], x0 + xs[max(idx)] + COL_W
        cv.line([a, b], [y + LANE_H / 2] * 2, color="#9aa1aa", lw=0.8)
        for ci in idx:
            cv.rect(x0 + xs[ci] + 0.5, y, COL_W - 1, LANE_H, fc="#3b4552", lw=0)
    return (xs[-1] + COL_W) if xs else COL_W


def _logo_height(al):
    return 12 + max(48, min(96, al["n"] * 8)) + 12 + 4


def _draw_partner_logo(cv, al, x0, y0):
    ph = max(48, min(96, al["n"] * 8))
    y_logo, y_off = y0 + 12, y0 + 12 + ph + 12
    cv.text(x0 - 8, y0 + 6, f"k of {al['n']}", size=6.5, color=INK3, ha="right", va="center")
    cv.text(x0 - 8, y_logo + ph / 2, "partner residues", size=6.5, color=INK3, ha="right", va="center")
    cv.text(x0 - 8, y_logo + ph + 7, "offset", size=6.5, color=INK3, ha="right", va="center")
    x = x0
    for c in al["columns"]:
        if c["kind"] == "gap":
            cv.text(x + GAP_W / 2, y_logo + ph - 4, "⋯", size=9, color=INK3, ha="center")
            cv.text(x + GAP_W / 2, y_off, f"{c['skipped']} aa", size=6, color=INK3, ha="center")
            x += GAP_W
            continue
        if c["offset"] == 0:
            cv.rect(x, y_logo - 1, PCOL_W, ph + 2, fc="#eef0f3", lw=0)
        cv.text(x + PCOL_W / 2, y0 + 9, str(c["k"]) if c["k"] else "·", size=6, ha="center",
                color=INK2 if c["k"] else "#c3c7cd")
        y = y_logo + ph
        for l, n in reversed(c["counts"]):
            h = n / max(1, al["n"]) * ph
            cv.glyph(l, x, y, PCOL_W, h, aa_color(l))
            y -= h
        cv.text(x + PCOL_W / 2, y_off, f"+{c['offset']}" if c["offset"] > 0 else str(c["offset"]), size=6,
                color=INK if c["offset"] == 0 else INK3, ha="center", weight="bold" if c["offset"] == 0 else "normal")
        x += PCOL_W
    return x - x0


def _wrap(cv, items, x0, y, max_w, size=7, gap=10):
    """Lay out text chips left to right, wrapping at x0 + max_w; items = [(text, color, weight)].
    Widths are measured on the rendered text (1 px = 1 unit at dpi 100). Returns the next y."""
    r = cv.fig.canvas.get_renderer()
    x = x0
    for t, c, w in items:
        h = cv.ax.text(x, y, t, fontsize=size, color=c, fontweight=w)
        tw = h.get_window_extent(renderer=r).width
        if x > x0 and x + tw > x0 + max_w:
            x, y = x0, y + 13
            h.set_position((x, y))
        x += tw + gap
    return y + 15


def interface_card(doc, iface_index, path, dpi=200):
    """One interface as the website's card shows it with every block open."""
    i = doc["interfaces"][iface_index]
    entry = i["motif"]
    name = {p["id"]: label(p["identity"], p["id"]) for p in i["partners"] or []}
    W_min = 760
    blocks_h = 120
    if entry:
        bait_w = GUTTER + 20 + len(entry["panel1"]) * (COL_W + JUMP_W)
        blocks_h += _bait_side_height(entry) + 40
        for m in entry["motifs"]:
            blocks_h += 40 + _logo_height(entry["alignments"][m["id"]]) + 16 * (2 + len(m["carriers"]) // 6) \
                + 14 * (len(m["elements"]) + 2) + 20
    else:
        bait_w = W_min
    widths = [bait_w] + ([GUTTER + 20 + sum(GAP_W if c["kind"] == "gap" else PCOL_W
                                             for c in entry["alignments"][m["id"]]["columns"])
                          for m in entry["motifs"]] if entry else [])
    W = int(max(W_min, max(widths) + 30))
    cv = _Canvas(W, int(blocks_h + 40))
    c = iface_color(iface_index)
    cv.rect(16, 12, 14, 14, fc=c if i["significant"] else c + "38", ec=c, lw=0 if i["significant"] else 0.8,
            hatch=None if i["significant"] else "////")
    n_p = len(i["partners"] or [])
    cv.text(36, 24, f"{_bait_name(doc)} {i['id']}", size=11, weight="bold")
    cv.text(16, 42, f"bait {compact(i['residues'])}  ·  {'significant' if i['significant'] else 'not significant'}"
            f"  ·  {n_p} partner{'s' if n_p != 1 else ''}", size=8, color=INK2)
    y = 66
    if not entry:
        cv.text(16, y, "No partner the motif engine could use touches this interface.", size=8, color=INK2)
        cv.save(path, dpi)
        return
    p1 = entry["panel1"]
    n_t = sum(1 for x in p1 if x["q"] is not None)
    n_s = sum(1 for x in p1 if x["q"] is not None and x["q"] < Q_SIG)
    cv.text(16, y, "Motifs — bait side", size=9, weight="bold")
    y += 15
    cv.text(16, y, f"{len(p1)} bait residues · {n_t} tested · {n_s} with q < {Q_SIG} · partners used "
            f"{entry['partners_used']}/{entry['partners_total']} · shuffle control: "
            f"{entry['control']['false_elements']} false element(s), columns "
            f"{'/'.join(map(str, entry['control']['false_columns']))}", size=7, color=INK2)
    y += 12
    if p1:
        _draw_bait_side(cv, entry, doc["bait"]["sequence"], GUTTER + 16, y)
    y += _bait_side_height(entry) + 8
    if not entry["motifs"]:
        cv.text(16, y + 10, "No motif found in this interface (no element passed the engine's tests) — "
                "a valid result, not missing data.", size=8, color=INK2)
    elif len(entry["motifs"]) > 1:
        rel = " · ".join(f"{r['a']}–{r['b']} {'same pocket' if r['relation'] == 'same_pocket' else 'separate sites'}"
                         f" · {'different (competing)' if r['status'] == 'different' else 'unresolved, too few carriers to tell'}"
                         for r in entry["motif_relations"])
        cv.text(16, y + 10, rel, size=7.5, color=INK2)
    y += 26
    for m in entry["motifs"]:
        al = entry["alignments"][m["id"]]
        cv.line([16, W - 16], [y - 8, y - 8], color=LINE, lw=0.8)
        head = (f"{len(m['positions'])} bait res · {m['n_groups']} group{'s' if m['n_groups'] != 1 else ''} · "
                f"{len(m['carriers'])} carrier{'s' if len(m['carriers']) != 1 else ''} · "
                f"{len(m['elements'])} element{'s' if len(m['elements']) != 1 else ''}"
                + (" · bipartite" if m["bipartite"] else "")
                + "".join(f" · linked by {x}" for x in m["linked_by"]))
        cv.text(16, y + 4, m["id"], size=9, weight="bold")
        cv.text(44, y + 4, head, size=7.5, color=INK2)
        note = f"Aligned on bait residue {al['anchor']}."
        if any(x["dir"] == -1 for x in al["carriers"]):
            note += " ← carriers bind antiparallel and are read backwards (C→N)."
        if al["orientation"] == "per-carrier":
            note += " Direction set per carrier from its own contacts."
        if al["n"] < len(al["carriers"]):
            note += f" {len(al['carriers']) - al['n']} carrier(s) share no bait position with the anchor holders (not in logo)."
        if al["missing_letters"]:
            note += f" {al['missing_letters']} residue(s) outside the listed ranges show as ?."
        cv.text(16, y + 18, note, size=7, color=INK3)
        y += 26
        _draw_partner_logo(cv, al, GUTTER + 16, y)
        y += _logo_height(al) + 6
        groups = {}
        for car in al["carriers"]:
            groups.setdefault(car["group"], []).append(car)
        chips = []
        for g in groups.values():
            t = "  ".join(f"{name.get(x['id'], x['id'])} {'→' if x['dir'] == 1 else '←'}"
                          + ("" if x["aligned"] else " (not in logo)") for x in g)
            chips.append((f"[{t}]" if len(g) > 1 else t, INK, "normal"))
        cv.text(16, y, "carriers", size=6.5, color=INK3)
        y = _wrap(cv, chips, 60, y, W - 80)
        cv.text(16, y, "el.     bait residues                         groups   q          or.   spacing   carriers",
                size=6.5, color=INK3, family="monospace")
        y += 13
        arrow = {"parallel": "→", "antiparallel": "←", "mixed": "↔"}
        for ei in m["elements"]:
            e = entry["elements"][ei]
            pos = compact(e["positions"])
            cv.text(16, y, f"E{ei + 1:<6}{pos[:36]:<38}{e['n_groups']:<9}{fmt_q(e['q']):<11}"
                    f"{arrow.get(e['orientation'], '?'):<6}{round((e['consistency'] or 0) * 100):>3}%      "
                    f"{len(e['carriers'])}", size=6.5, color=INK, family="monospace")
            y += 13
        y += 22
    cv.text(16, y, f"Motif engine {doc['params'].get('motif_engine')} · contact Cβ–Cβ ≤ "
            f"{doc['params'].get('motif_touch_cb')} Å · partners ≥ {round(100 * doc['params'].get('motif_id_group', 0))}% "
            f"identical = one group · {doc['params'].get('motif_n_perm')} letter shuffles", size=6.5, color=INK3)
    cv.text(16, y + 12, f"column height ∝ √ groups (busiest = full) · ✱ q < {Q_SIG} vs letter shuffles · hatched = "
            "untested (< 3 groups) · → carrier read N→C, ← read backwards (antiparallel)", size=6.5, color=INK3)
    cv.crop(y + 24)
    cv.save(path, dpi)


def partner_strips(doc, iface_index, path, dpi=200, tiers=("H", "M", "L")):
    """Every prey at the interface (rows in the annotation's order: most H residues first)
    with its residue ranges as tier-styled chips; residues of any motif it carries in bold."""
    i = doc["interfaces"][iface_index]
    parts = i["partners"] or []
    motif_res = {}
    for m in ((i["motif"] or {}).get("motifs") or []):
        for car in i["motif"]["alignments"][m["id"]]["carriers"]:
            motif_res.setdefault(car["id"], {})[m["id"]] = set(car["residues"])
    W, ROW, X0 = 1000, 18, 205

    def layout(a):
        """Chip positions for one partner, wrapping onto extra lines: [(line, x, w, text, tier, bold)]."""
        mres = set().union(*motif_res[a["id"]].values()) if a["id"] in motif_res else set()
        out, line, x = [], 0, X0
        for s_, e, t, seq in a["ranges"]:
            if t not in tiers:
                continue
            txt = f"{s_}–{e} {seq}" if e > s_ else f"{s_} {seq}"
            w = 5.6 * len(txt) + 10
            if x > X0 and x + w > W - 10:
                line, x = line + 1, X0
            out.append((line, x, w, txt, t, any(r in mres for r in range(s_, e + 1))))
            x += w + 4
        return out
    lays = [layout(a) for a in parts]
    H = 70 + sum(ROW * (1 + (max((c[0] for c in L), default=0))) for L in lays) + 20
    cv = _Canvas(W, H)
    cv.text(16, 22, f"{_bait_name(doc)} {i['id']} — partners ({len(parts)})", size=10, weight="bold")
    cv.text(16, 38, "every partner with ≥ 1 contact at PAE ≤ 10 Å; its residues tiered by their best contact: "
            "H PAE ≤ 5 Å (precision 72 %), M 5–10 Å (37 %), L > 10 Å (14 %); bold = motif residues",
            size=7, color=INK3)
    cv.text(150, 54, "H/M/L", size=6.5, color=INK3, family="monospace")
    if i["partners"] is None:
        cv.text(16, 70, "No pairwise PAE for this screen: the partner side is unknown (not empty).", size=8)
    y = 58
    for a, L in zip(parts, lays):
        tag = (" " + ",".join(sorted(motif_res[a["id"]]))) if a["id"] in motif_res else ""
        cv.text(16, y + 11, label(a["identity"], a["id"]) + tag, size=7, weight="bold")
        cv.text(150, y + 11, f"{a['n_high']}/{a['n_medium']}/{a['n_low']}", size=6.5, color=INK3, family="monospace")
        for line, x, w, txt, t, bold in L:
            st = TIER_STYLE[t]
            yy = y + line * ROW
            cv.rect(x, yy + 2, w, 13, fc=st["fc"], ec=st["ec"], ls=st["ls"], lw=0.7)
            cv.text(x + 4, yy + 11.5, txt, size=6.5, color=st["tc"], family="monospace",
                    weight="bold" if bold else "normal")
        y += ROW * (1 + max((c[0] for c in L), default=0))
    cv.save(path, dpi)
