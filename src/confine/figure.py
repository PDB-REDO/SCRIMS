"""
The matrix figure: partners x bait residues, with the called sites marked.

Design decisions worth knowing before changing them:

* **Sequential, one hue, light -> dark** for the signal (magnitude). No rainbow:
  a rainbow ramp invents category boundaries where the data has none.
* **0 is not a pale value, it is a different statement.** "No contact in this
  prediction" is categorical, not a small magnitude, so zeros are masked and
  drawn in a neutral grey with their own legend entry. Colouring them as the
  bottom of the ramp would make absence look like weak evidence.
* **The colourbar is labelled in PAE, not in signal units**, because the signal
  is a transform nobody reads natively. 0.5 is exactly PAE = d0_signal (10 A).
* **Sites are marked by position and label, never by colour alone** — a span
  bracket plus its residue range, with recurrent and strength-only spans given
  different line styles as well as different colours.

matplotlib is imported inside the functions, so it stays an optional extra.
"""

from __future__ import annotations

import numpy as np

NO_CONTACT = "#eceef0"
RAMP = "Blues"
INK = "#1a1a1a"
MUTED = "#6b7280"
REC = "#08519c"
STR = "#a63603"


def plot_matrix(res, ax=None, order=None, max_rows=None, show_sites=True,
                xlabel=True):
    """Draw one screen's matrix. `res` is a confine.Result.

    Returns the matplotlib Axes. With `max_rows`, the strongest rows by total
    signal are kept and the caption count is adjusted — stated on the figure, so
    a truncated matrix can never be mistaken for a complete one.
    """
    import matplotlib.pyplot as plt

    M = res.M_signal
    if M.shape[0] == 0:
        raise ValueError("no rows in the matrix — nothing to plot")
    if order is None:
        order, _ = res.row_order()
    rows = list(order)
    truncated = 0
    if max_rows and len(rows) > max_rows:
        keep = set(np.argsort(M.sum(1))[::-1][:max_rows].tolist())
        truncated = len(rows) - len(keep)
        rows = [i for i in rows if i in keep]

    Mv = np.ma.masked_where(M[rows] <= 0, M[rows])
    if ax is None:
        _fig, ax = plt.subplots(figsize=(9, max(2.2, 0.06 * len(rows) + 1.4)))
    # `with_extremes` rather than `set_bad`: the latter is pending deprecation
    # in matplotlib 3.11 and this package supports >= 3.7, where both exist.
    cmap = plt.get_cmap(RAMP).with_extremes(bad=NO_CONTACT)
    im = ax.imshow(Mv, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                   interpolation="nearest",
                   extent=[res.res_ids[0] - 0.5, res.res_ids[-1] + 0.5, len(rows) - 0.5, -0.5])
    if xlabel:  # in the stacked figure the track below carries the axis
        ax.set_xlabel(f"{res.bait.get('gene', 'bait')} residue")
    ax.set_ylabel(f"partner predictions (n = {len(rows)})"
                  + (f", {truncated} not shown" if truncated else ""))
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelcolor=INK)

    if show_sites:
        _mark_sites(ax, res, len(rows))
    return ax, im


def _mark_sites(ax, res, n_rows):
    """Site spans as brackets under the matrix, labelled with their residue range."""
    y = n_rows - 0.5
    for s in res.sites:
        if not s["significant"]:
            continue
        rec = bool(s.get("recurrent"))
        ax.plot([s["res_start"], s["res_end"]], [y + 0.6, y + 0.6],
                lw=2.4, solid_capstyle="butt",
                color=REC if rec else STR, linestyle="-" if rec else (0, (2, 1.2)),
                clip_on=False, zorder=5)


def figure(res, path=None, max_rows=None, dpi=300, height=None, width=9.0):
    """The package's headline output: matrix + summed-signal track + site table row.

    Writes to `path` if given (extension picks the format: .png, .pdf, .svg) and
    returns the Figure.

    Height defaults to one that grows with the number of partners; set `height`
    (inches) for a fixed size, e.g. to fit a journal column. Rows simply get
    thinner — no partner is dropped, unlike `max_rows`.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    n = res.M_signal.shape[0] if max_rows is None else min(max_rows, res.M_signal.shape[0])
    h = height if height else max(3.4, 0.055 * n + 2.6)
    fig, (ax, axt) = plt.subplots(
        2, 1, figsize=(width, h), sharex=True,
        gridspec_kw=dict(height_ratios=[max(1.4, h - 1.4), 1.0], hspace=0.08))

    _ax, im = plot_matrix(res, ax=ax, max_rows=max_rows, xlabel=False)
    cb = fig.colorbar(im, ax=[ax, axt], fraction=0.025, pad=0.015)
    cb.set_label("confidence signal   1/(1+(PAE/10)²)")
    cb.set_ticks([0.0, 0.2, 0.5, 0.8, 1.0])
    cb.ax.set_yticklabels(["0", "0.2", "0.5\n(PAE 10 Å)", "0.8\n(PAE 5 Å)", "1"])
    cb.outline.set_visible(False)

    # Track: total signal per residue — what the recurrence statistic tests.
    axt.fill_between(res.res_ids, res.s_sum, color=REC, alpha=0.18, linewidth=0)
    axt.plot(res.res_ids, res.s_sum, color=REC, lw=1.4)
    axt.set_ylabel("summed\nsignal")
    axt.set_xlabel(f"{res.bait.get('gene', 'bait')} residue  (length {res.bait['length']})")
    axt.margins(x=0)
    for s in ("top", "right"):
        axt.spines[s].set_visible(False)
    for s in ("bottom", "left"):
        axt.spines[s].set_color(MUTED)
    axt.tick_params(colors=MUTED, labelcolor=INK)

    for s in res.sites:
        if not s["significant"]:
            continue
        rec = bool(s.get("recurrent"))
        axt.axvspan(s["res_start"] - 0.5, s["res_end"] + 0.5, color=REC if rec else STR,
                    alpha=0.13, linewidth=0, zorder=0)

    handles = [
        Line2D([], [], color=REC, lw=2.4, label="recurrent site (q < 0.05)"),
        Line2D([], [], color=STR, lw=2.4, linestyle=(0, (2, 1.2)),
               label="strength site (criterion, no p-value)"),
        Line2D([], [], color=NO_CONTACT, lw=6, label="no contact (0)"),
    ]
    axt.legend(handles=handles, loc="upper left", bbox_to_anchor=(0, -0.42), ncol=3,
               frameon=False, fontsize=8, labelcolor=INK)

    if path:
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return fig
