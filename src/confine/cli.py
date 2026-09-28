"""
Command line: one screen in, a figure and the called sites out.

    confine run SCREEN_DIR -o OUTDIR

Writes, in OUTDIR:
    matrix.png       the figure: matrix + summed-signal track + marked sites
    sites.tsv        one row per called site, with significance, channel and interface
    interfaces.tsv   one row per 3-D interface, with its member segments
    matrix.npz       M_pae (raw Angstrom, NaN = no contact), M_signal, res_ids, preds
    summary.json     bait, parameters, counts, and every skipped prediction with a reason
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from confine import DEFAULTS, analyze

SITE_COLUMNS = ["res_start", "res_end", "interface", "channel", "significant",
                "recurrent", "q_recur", "recurrent_core", "strength_significant",
                "strength_core", "n_partners_strong", "strength_A", "strength_A_per_res"]


def _write_sites(res, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(SITE_COLUMNS)
        for s in res.sites:
            w.writerow([_fmt(s.get(c)) for c in SITE_COLUMNS])


def _write_interfaces(res, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["interface", "n_patches", "significant", "channels", "n_residues", "segments"])
        for _, r in res.interfaces.iterrows():
            w.writerow([r.interface, r.n_patches, bool(r.significant), r.channels,
                        len(r.residues),
                        ";".join(f"{a}-{b}" for a, b in r.segments)])


def _fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    if isinstance(v, (list, tuple)):
        return ";".join(f"{a}-{b}" for a, b in v) if v and isinstance(v[0], (list, tuple)) else \
               ";".join(str(x) for x in v)
    return v


def main(argv=None):
    ap = argparse.ArgumentParser(prog="confine", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    run = sub.add_parser("run", help="analyse one screen directory")
    run.add_argument("screen_dir", type=Path)
    run.add_argument("-o", "--out", type=Path, default=Path("confine_out"))
    run.add_argument("--no-figure", action="store_true", help="skip the figure")
    run.add_argument("--figure-format", default="png", choices=["png", "pdf", "svg"])
    run.add_argument("--max-rows", type=int, default=None,
                     help="cap the plotted rows (the figure states how many are hidden)")
    run.add_argument("--bait-chain", dest="bait_chain", default=None,
                     help="name the bait chain instead of detecting it by sequence identity")
    for k, v in DEFAULTS.items():
        if k == "bait_chain":
            continue
        if isinstance(v, bool):
            run.add_argument(f"--{k.replace('_', '-')}", dest=k,
                             action="store_false" if v else "store_true", default=v)
        elif v is None:
            run.add_argument(f"--{k.replace('_', '-')}", dest=k, type=float, default=None)
        else:
            run.add_argument(f"--{k.replace('_', '-')}", dest=k, type=type(v), default=v)
    a = ap.parse_args(argv)

    params = {k: getattr(a, k) for k in DEFAULTS}
    res = analyze(a.screen_dir, **params)
    a.out.mkdir(parents=True, exist_ok=True)

    _write_sites(res, a.out / "sites.tsv")
    _write_interfaces(res, a.out / "interfaces.tsv")
    np.savez_compressed(a.out / "matrix.npz", M_pae=res.M_pae, M_signal=res.M_signal,
                        res_ids=res.res_ids, preds=np.array(res.preds, dtype=object),
                        allow_pickle=True)
    n_sig = len(res.significant())
    (a.out / "summary.json").write_text(json.dumps({
        "screen_dir": str(a.screen_dir),
        "bait": {k: v for k, v in res.bait.items() if k != "chain_ids"},
        "n_predictions_in_matrix": len(res.preds),
        "n_sites": len(res.sites), "n_significant": n_sig,
        "n_recurrent": sum(1 for s in res.sites if s.get("recurrent")),
        "n_strength": sum(1 for s in res.sites if s.get("strength_significant")),
        "n_interfaces": len(res.interfaces),
        "single_pair": res.single_pair,
        "skipped": res.skipped, "params": params,
    }, indent=2, default=str))

    if not a.no_figure and res.M_signal.shape[0]:
        from confine.figure import figure
        fig = figure(res, a.out / f"matrix.{a.figure_format}", max_rows=a.max_rows)
        try:
            import matplotlib.pyplot as plt
            plt.close(fig)
        except Exception:
            pass

    print(f"{len(res.preds)} predictions -> {len(res.sites)} sites "
          f"({n_sig} significant) in {len(res.interfaces)} interfaces; wrote {a.out}")
    if res.skipped:
        print(f"skipped {len(res.skipped)} predictions (see summary.json)")
    if res.single_pair:
        print("single prediction: recurrence not tested, strength regions only")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
