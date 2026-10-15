"""
Static figures of the companion website's views, drawn from an `scrims/1`
document (scrims.export.document), so the figures, the JSON and the HTML report
show the same thing.

  overview(doc, path)            confidence matrix (binder rows, clustered; contacted
                                 residues as columns) + RMSF and recurrence tracks +
                                 one interface lane
  interface_card(doc, i, path)   one interface: bait-side view in the website's default
                                 "evidence view" (Panel 1 letter columns, height ~
                                 sqrt(groups), * q < 0.05, hatched = untested, the bait
                                 residue as a Zappo box), motif lanes, the relations line,
                                 then per motif its partner-sequence logo in the default
                                 "classic" view (bits, group-weighted, exact small-sample
                                 correction, scaled to its tallest column) over the
                                 per-carrier table in the same columns
  partner_strips(doc, i, path)   every prey at the interface with its tiered residue
                                 ranges (H solid, M light, L dashed), motif residues bold

Colours, sizes and rules follow the website (webiste_viz_int src/lib/{colors,logo}.ts and
src/components/{Heatmap,MotifEngineView}.tsx, 2026-10-07): Greens ramp for confidence,
purple ramp for RMSF, eight categorical interface hues keyed by interface order
(significant = solid, not = hatched), amino acids in the official Jalview Zappo colours
(logo letters filled + thin dark outline; small letters a coloured box + black letter),
selection highlight cyan. Logo letters are stacked contiguously from the baseline,
largest on top, each fitted to its slot from its exact outline; shares too small to
draw form one grey band on top (scrims.logo). Colour is never the only cue: every
band, column and chip also carries text.

Layout is in pixels on one axes per figure (1 px = 1 data unit at dpi 100); `dpi`
sets the output resolution.
"""

from __future__ import annotations

import math

import numpy as np

from scrims.export import label
from scrims.logo import bits_axis_max, build_stack, group_weighted, info_content

CONF_LO, CONF_HI = "#f7fcf5", "#00441b"
RMSF_LO, RMSF_HI = "#efedf5", "#3f007d"
def motif_label(mid):
    """'M1' -> 'Motif 1' (2026-10-12, user's PI: "M1 is not obvious"; as the report and website); data ids stay M1."""
    return f"Motif {mid[1:]}" if mid[:1] == "M" and mid[1:].isdigit() else mid


HIGHLIGHT = "#00ffff"          # selection colour, as the website's (gold collided with Zappo aromatic)
INTERFACE_COLORS = ["#2a78d6", "#d06f00", "#d4609a", "#a0662f", "#1a9fc0", "#d4462c", "#7b61d9", "#c2566b"]
# official Jalview Zappo (jalview.org/help/html/colourSchemes/zappo.html)
AA_GROUPS = [("aliphatic / hydrophobic", "ILVAM", "#ffafaf"), ("aromatic", "FWY", "#ffc800"),
             ("positive", "KRH", "#6464ff"), ("negative", "DE", "#ff0000"),
             ("hydrophilic", "STNQ", "#00ff00"), ("conformationally special", "PG", "#ff00ff"),
             ("cysteine", "C", "#ffff00")]
LETTER_OUTLINE = "#1f2329"
OTHER_FILL, GREY_LETTER = "#c3c7cd", "#d3d7dc"
AA_COLOR = {a: c for _, letters, c in AA_GROUPS for a in letters}
INK, INK2, INK3, LINE = "#191c20", "#4a5058", "#7c828b", "#e2e4e8"
Q_SIG = 0.05
TIER_STYLE = {"H": dict(fc="#25603c", ec="#25603c", ls="-", tc="white"),     # the website's tier chips
              "M": dict(fc="#adc5b4", ec="#adc5b4", ls="-", tc=INK),
              "L": dict(fc="white", ec=INK3, ls="--", tc=INK2)}


def aa_color(ch):
    return AA_COLOR.get(ch.upper(), "#c3c7cd")


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
        """Letter `ch` with its exact outline fitted into the w x h box whose bottom edge is
        at y_bottom (the website's Glyph: height exact; width 86 % of the column for an
        H-wide letter, wide letters squeezed to at most 94 %), filled + thin dark outline."""
        from matplotlib.font_manager import FontProperties
        from matplotlib.patches import PathPatch
        from matplotlib.textpath import TextPath
        from matplotlib.transforms import Affine2D
        fp = FontProperties(family="DejaVu Sans", weight="bold")
        p = TextPath((0, 0), ch, size=10, prop=fp)
        cap = TextPath((0, 0), "H", size=10, prop=fp).get_extents()
        bb = p.get_extents()
        sx = min(0.86 * w / max(cap.width, 1e-6), 0.94 * w / max(bb.width, 1e-6))
        sy = h / max(bb.height, 1e-6)
        t = (Affine2D().translate(-bb.x0 - bb.width / 2, -bb.y0).scale(sx, -sy)
             .translate(x + w / 2, y_bottom))
        self.ax.add_patch(PathPatch(t.transform_path(p), fc=color, ec=LETTER_OUTLINE,
                                    lw=0.25 if h < 8 else 0.4))

    def stack(self, letters, x, y_base, w, grey=False):
        """[(letter, height)] summing to the column height, baseline at y_base: visible letters
        contiguous from the baseline, largest on top; too-small shares pooled into one grey band."""
        vis, other, _ = build_stack(letters)
        y = y_base
        for l, h in vis:
            self.glyph(l, x, y, w, h, GREY_LETTER if grey else aa_color(l))
            y -= h
        if other > 0.3:
            self.rect(x + 2, y - other, w - 4, other, fc=OTHER_FILL, lw=0)

    def aabox(self, ch, x, y_top, w, h, size=7):
        """A small residue, Jalview-style: Zappo box, black letter."""
        self.rect(x, y_top, w, h, fc=aa_color(ch), lw=0)
        self.text(x + w / 2, y_top + h / 2, ch, size=size, color="#111418", ha="center", va="center",
                  weight="bold", family="monospace")

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
PCOL_W, GAP_W, GUTTER, PLEFT, PLOGO_H, CROW_H = 18, 30, 110, 210, 80, 16


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
        tot = sum(c["letters"].values()) or 1
        cv.stack([(l, v / tot * col_h) for l, v in c["letters"].items()], X, y_logo + LOGO_H, COL_W, grey=untested)
        b = seq[c["pos"] - 1] if seq and c["pos"] - 1 < len(seq) else "?"
        cv.aabox(b, X + 1.5, y_bait + 2, COL_W - 3, 12)
        cv.text(X + COL_W / 2, y_num + 2, str(c["pos"]), size=6, color=INK3, ha="center", va="top",
                rotation=90, family="monospace")
    col_of = {c["pos"]: k for k, c in enumerate(cols)}
    for li, m in enumerate(entry["motifs"]):
        y = y_lanes + li * (LANE_H + LANE_GAP)
        idx = [col_of[p] for p in m["positions"] if p in col_of]
        cv.text(x0 - 8, y + LANE_H / 2, motif_label(m["id"]), size=7, color=INK, ha="right", va="center", weight="bold")
        if not idx:
            continue
        a, b = x0 + xs[min(idx)], x0 + xs[max(idx)] + COL_W
        cv.line([a, b], [y + LANE_H / 2] * 2, color="#9aa1aa", lw=0.8)
        for ci in idx:
            cv.rect(x0 + xs[ci] + 0.5, y, COL_W - 1, LANE_H, fc="#3b4552", lw=0)
    return (xs[-1] + COL_W) if xs else COL_W


def _partner_cols(al):
    """Per logo column: ([(letter, f)], n carriers, bits), or None for a gap column. Every carrier counts
    (viz_int 2026-10-08: paralogs counted once hid VAPA FFAT's second F), as on the website."""
    out = []
    for c in al["columns"]:
        if c["kind"] == "gap":
            out.append(None)
            continue
        o = str(c["offset"])
        items = [(car["cells"][o], car["id"]) for car in al["carriers"] if o in (car.get("cells") or {})]
        freq, n = group_weighted(items)
        out.append((freq, n, info_content([f for _, f in freq], n)))
    return out


def _logo_height(al):
    """Height _draw_partner_logo uses (axis note + groups row + logo + offsets); kept for callers
    that lay out the logo themselves (paper/make_figure1.py)."""
    return 16 + 12 + PLOGO_H + 12 + 6


def _logo_width(al):
    return sum(GAP_W if c["kind"] == "gap" else PCOL_W for c in al["columns"])


def _draw_partner_logo(cv, al, x0, y0):
    """The website's default partner-side view: classic logo in bits, y axis scaled to the
    tallest column (rounded up to 0.5 bit). Returns the height used."""
    pc = _partner_cols(al)
    ics = [d[2] if d else 0.0 for d in pc]
    amax = bits_axis_max(ics)
    empty = any(d is not None for d in pc) and all(d is None or d[1] <= 1 for d in pc)
    cv.text(16, y0 + 8, "classic logo is empty: every column comes from a single carrier "
            "(one sequence carries no information)" if empty else
            f"classic logo scaled to the tallest column (max {max(ics + [0]):.2f} bits; axis 0–{amax:.1f} bits)"
            + ("; antiparallel to the bait, shown N→C" if al.get("flipped") else ""),
            size=6.5, color=INK3)
    y0 += 16
    y_logo, y_off = y0 + 12, y0 + 12 + PLOGO_H + 12
    cv.text(x0 - 8, y0 + 6, "carriers", size=6.5, color=INK3, ha="right", va="center")
    cv.text(x0 - 34, y_logo + PLOGO_H / 2, "bits", size=6.5, color=INK3, ha="right", va="center")
    step = 0.5 if amax <= 2.5 else 1.0
    b = 0.0
    while b <= amax + 1e-9:
        cv.text(x0 - 8, max(y_logo + 3, y_logo + PLOGO_H - b / amax * PLOGO_H), f"{b:.1f}" if step < 1 else f"{b:g}",
                size=5.5, color=INK3, ha="right", va="center", family="monospace")
        b += step
    cv.text(x0 - 8, y_logo + PLOGO_H + 7, "offset", size=6.5, color=INK3, ha="right", va="center")
    cv.line([x0, x0 + _logo_width(al)], [y_logo + PLOGO_H] * 2, color="#d0d3d9", lw=0.8)
    x = x0
    for c, d in zip(al["columns"], pc):
        if c["kind"] == "gap":
            cv.text(x + GAP_W / 2, y_logo + PLOGO_H - 4, "⋯", size=9, color=INK3, ha="center")
            cv.text(x + GAP_W / 2, y_off, f"{c['skipped']} aa", size=6, color=INK3, ha="center")
            x += GAP_W
            continue
        freq, n, ic = d
        cv.text(x + PCOL_W / 2, y0 + 9, str(n) if n else "·", size=6, ha="center", color=INK2 if n else "#c3c7cd")
        col_h = ic / amax * PLOGO_H
        cv.stack([(l, f * col_h) for l, f in freq], x, y_logo + PLOGO_H, PCOL_W)
        cv.text(x + PCOL_W / 2, y_off, f"+{c['offset']}" if c["offset"] > 0 else str(c["offset"]), size=6,
                color=INK if c["offset"] == 0 else INK3, ha="center", weight="bold" if c["offset"] == 0 else "normal")
        x += PCOL_W
    return _logo_height(al)


def _letters_of(partners):
    """prediction -> {partner residue: letter}, from the interface's partner ranges."""
    out = {}
    for p in partners or []:
        m = {}
        for a, b, _t, seq in p["ranges"]:
            for r in range(a, b + 1):
                m[r] = seq[r - a] if r - a < len(seq) else "?"
        out[p["id"]] = m
    return out


def _draw_carrier_table(cv, al, x0, y0, name, uni, seqs):
    """One row per carrier, in the logo's columns: gene, UniProt, paralog group (gN when the
    group has > 1 carrier here), reading arrow, residue range, its letter at each offset
    ('·' when absent) and its actual segment with the motif residues bold. Returns height."""
    gno, gsz = {}, {}
    for c in al["carriers"]:
        gno.setdefault(c["group"], len(gno) + 1)
        gsz[c["group"]] = gsz.get(c["group"], 0) + 1
    rows = sorted(al["carriers"], key=lambda c: gno[c["group"]])
    seg_x = x0 + _logo_width(al) + 10
    cv.text(16, y0 + 9, f"Carriers ({len(rows)})", size=7, color=INK2, weight="bold")
    y = y0 + 14
    for car in rows:
        cv.line([16, seg_x + 6 * (car["residues"][-1] - car["residues"][0] + 1) + 10], [y, y], color=LINE, lw=0.5)
        yc = y + CROW_H / 2
        nm = name.get(car["id"], car["id"])
        cv.text(16, yc, nm[:13], size=7, color=INK, va="center", weight="bold")
        u = uni.get(car["id"])
        if u and u != nm:
            cv.text(84, yc, u, size=5.8, color=INK3, va="center", family="monospace")
        if gsz[car["group"]] > 1:
            cv.text(126, yc, f"g{gno[car['group']]}", size=5.8, color=INK2, va="center", family="monospace",
                    bbox=dict(boxstyle="round,pad=0.15", fc="none", ec=INK3, ls="--", lw=0.5))
        cv.text(146, yc, "→" if car["dir"] == 1 else "←", size=7.5, color=INK2, va="center", weight="bold")
        lo, hi = car["residues"][0], car["residues"][-1]
        cv.text(x0 - 6, yc, str(lo) if lo == hi else f"{lo}–{hi}", size=6, color=INK2, ha="right", va="center",
                family="monospace")
        x = x0
        cells = car.get("cells") or {}
        for c in al["columns"]:
            if c["kind"] == "gap":
                x += GAP_W
                continue
            l = cells.get(str(c["offset"]))
            if l:
                cv.aabox(l, x + 2, y + 2, PCOL_W - 4, CROW_H - 4, size=6)
            else:
                cv.text(x + PCOL_W / 2, yc, "·", size=7, color="#c3c7cd", ha="center", va="center")
            x += PCOL_W
        mres, sq = set(car["residues"]), seqs.get(car["id"], {})
        for k, r in enumerate(range(lo, hi + 1)):
            l = sq.get(r)
            cv.text(seg_x + 6 * k, yc, "·" if l is None else l, size=6.5, va="center", family="monospace",
                    color="#c3c7cd" if l is None else ("#111418" if r in mres else INK3),
                    weight="bold" if r in mres else "normal")
        y += CROW_H
    return y - y0 + 4


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
    uni = {p["id"]: (p["identity"] or {}).get("uniprot") for p in i["partners"] or []}
    seqs = _letters_of(i["partners"])
    W_min = 760
    blocks_h = 120
    if entry:
        bait_w = GUTTER + 20 + len(entry["panel1"]) * (COL_W + JUMP_W)
        blocks_h += _bait_side_height(entry) + 40
        for m in entry["motifs"]:
            blocks_h += 40 + 16 + 12 + PLOGO_H + 30 + 20 + CROW_H * len(entry["alignments"][m["id"]]["carriers"]) + 30
    else:
        bait_w = W_min

    def _seg_len(al):
        return max([c["residues"][-1] - c["residues"][0] + 1 for c in al["carriers"]] or [0])
    widths = [bait_w] + ([16 + PLEFT + _logo_width(entry["alignments"][m["id"]]) + 20
                          + 6 * _seg_len(entry["alignments"][m["id"]]) for m in entry["motifs"]] if entry else [])
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
    cv.text(16, y, "Interface – Bait-side view", size=9, weight="bold")
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
        rel = " · ".join(f"{motif_label(r['a'])}–{motif_label(r['b'])} {'same pocket' if r['relation'] == 'same_pocket' else 'separate sites'}"
                         f" · {'different (competing)' if r['status'] == 'different' else 'unresolved, too few carriers to tell'}"
                         for r in entry["motif_relations"])
        cv.text(16, y + 10, rel, size=7.5, color=INK2)
    y += 26
    for m in entry["motifs"]:
        al = entry["alignments"][m["id"]]
        cv.line([16, W - 16], [y - 8, y - 8], color=LINE, lw=0.8)
        head = (f"{len(m['positions'])} bait res · {m['n_groups']} group{'s' if m['n_groups'] != 1 else ''} · "
                f"{len(m['carriers'])} carrier{'s' if len(m['carriers']) != 1 else ''}"
                + (" · bipartite" if m["bipartite"] else "")
                + "".join(f" · linked by {x}" for x in m["linked_by"]))
        cv.text(16, y + 4, motif_label(m["id"]), size=9, weight="bold")
        cv.text(66, y + 4, head, size=7.5, color=INK2)
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
        y += _draw_partner_logo(cv, al, 16 + PLEFT, y)
        y += _draw_carrier_table(cv, al, 16 + PLEFT, y, name, uni, seqs) + 16
    cv.text(16, y, f"Motif engine {doc['params'].get('motif_engine')} · contact Cβ–Cβ ≤ "
            f"{doc['params'].get('motif_touch_cb')} Å · partners ≥ {round(100 * doc['params'].get('motif_id_group', 0))}% "
            f"identical = one group · {doc['params'].get('motif_n_perm')} letter shuffles", size=6.5, color=INK3)
    for k, t in enumerate([
            "bait side (evidence view): column height ∝ √ groups (busiest = full), letters = shares",
            f"✱ q < {Q_SIG} vs letter shuffles · hatched = untested (< 3 groups) · grey band on top = letters too small to draw",
            "partner side (classic logo): bits = log2 20 − (H + e(n)), exact small-sample correction e(n); "
            "related partners count once",
            "→ carrier read N→C, ← read backwards (antiparallel, flipped) · bold = motif residues, grey · = outside "
            "the partner's listed ranges"]):
        cv.text(16, y + 12 * (k + 1), t, size=6.5, color=INK3)
    y += 48
    cv.text(16, y + 22, "Amino-acid colours:", size=6.5, color=INK2, va="center")
    x0_ = 178
    x, yl = x0_, y + 22
    for nm_, letters, _c in AA_GROUPS:
        need = 11 * len(letters) + 4 + 4.4 * len(nm_) + 12
        if x > x0_ and x + need > W - 16:
            x, yl = x0_, yl + 14
        for k, ch in enumerate(letters):
            cv.aabox(ch, x + 11 * k, yl - 7, 10, 10, size=5.5)
        x += 11 * len(letters) + 4
        cv.text(x, yl, nm_, size=6, color=INK3, va="center")
        x += 4.4 * len(nm_) + 12
    cv.crop(yl + 10)
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
