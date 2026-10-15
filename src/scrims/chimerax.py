"""
A ChimeraX command script per screen: the 3-D counterpart of the figures.

It opens the model the website's viewer shows by default (the first binder in the
matrix's display order), colours the bait by interface (the figures' interface
hues; grey elsewhere, partner chain light grey), names every interface and every
motif as a ChimeraX selection (`I3`, `I3_M1`), shows motif bait residues as sticks
and, per motif, opens its first carrier's model with that carrier's motif residues
named too (`I3_M1_<partner>`), so `view I3_M1` / `show I3_M1_PXN atoms` work at once.

Models kept only as `.cif.zst` are written next to the script (ChimeraX cannot
read zstd), under `models/`.
"""

from __future__ import annotations

import shutil
from pathlib import Path

from scrims.export import label
from scrims.views import HIGHLIGHT, iface_color


def _spec(model, chain, residues):
    rs = sorted(set(residues))
    out, i = [], 0
    while i < len(rs):
        j = i
        while j + 1 < len(rs) and rs[j + 1] == rs[j] + 1:
            j += 1
        out.append(f"{rs[i]}-{rs[j]}" if j > i else str(rs[i]))
        i = j + 1
    return f"#{model}/{chain}:{','.join(out)}" if out else None


def write_cxc(doc, model_path_of, outdir):
    """Write `<outdir>/scrims.cxc`. `model_path_of(pred)` returns a readable model CIF
    path (a plain file; temporary paths are copied into `<outdir>/models/`)."""
    outdir = Path(outdir)
    (outdir / "models").mkdir(exist_ok=True)

    def local(pred):
        src = Path(model_path_of(pred))
        dst = outdir / "models" / f"{pred}.cif"
        if not dst.exists():
            shutil.copyfile(src, dst)
        return f"models/{pred}.cif"

    rows = doc["matrix"]["row_order"]
    if not rows:
        return None
    rep = doc["predictions"][rows[0]]["id"]
    bc = doc["bait"]["chain_ids"][rep]
    bait = (doc["bait"].get("identity") or {}).get("gene") or doc["screen"]["name"]
    L = [f"# SCRIMS {doc['scrims_version']} -- {bait}: interfaces and motifs in 3-D",
         f"# model #1 = {rep} ({label(doc['predictions'][rows[0]]['identity'], rep)}), the first binder in display order",
         "# run from this folder:  chimerax scrims.cxc   (or `open scrims.cxc` inside ChimeraX)",
         f"open {local(rep)}",
         "set bgColor white", "color #1 #b9bfc7", "cartoon #1", "hide #1 atoms",
         f"transparency #1 & ~#1/{bc} 75 target c",      # the partner chain fades; the bait carries the colour
         "lighting soft", "graphics silhouettes true"]
    model = 1
    for idx, i in enumerate(doc["interfaces"]):
        sp = _spec(1, bc, i["residues"])
        if not sp:
            continue
        L += [f"# {i['id']}: bait {i['residues'][0]}..{i['residues'][-1]}, "
              f"{'significant' if i['significant'] else 'not significant'}, {len(i['partners'] or [])} partners",
              f"name {i['id']} {sp}", f"color {i['id']} {iface_color(idx)}"]
        for m in (i["motif"] or {}).get("motifs", []):
            msp = _spec(1, bc, m["positions"])
            L += [f"name {i['id']}_{m['id']} {msp}", f"show {i['id']}_{m['id']} atoms", f"style {i['id']}_{m['id']} stick",
                  f"color {i['id']}_{m['id']} {HIGHLIGHT} targ a"]
            al = i["motif"]["alignments"][m["id"]]
            car = next((c for c in al["carriers"] if c["aligned"]), None)
            if car is None:
                continue
            model += 1
            cbc = doc["bait"]["chain_ids"][car["id"]]
            ident = next((p["identity"] for p in i["partners"] or [] if p["id"] == car["id"]), None)
            nm = label(ident, car["id"]).replace(" ", "_")
            L += [f"# {i['id']} {m['id']} carrier {nm} ({car['id']}), motif residues {car['residues']}",
                  f"open {local(car['id'])}", f"matchmaker #{model}/{cbc} to #1/{bc}",
                  f"color #{model} #d9dce1", f"hide #{model} models"]
            # the partner chain is the model's non-bait chain; ChimeraX selects it as "not the bait chain"
            psp = f"#{model}:" + ",".join(map(str, car["residues"])) + f" & ~#{model}/{cbc}"
            L += [f"name {i['id']}_{m['id']}_{nm} {psp}", f"color {i['id']}_{m['id']}_{nm} {HIGHLIGHT}"]
    L += ["view #1", "# named selections: " + " ".join(
        [i["id"] for i in doc["interfaces"]] +
        [f"{i['id']}_{m['id']}" for i in doc["interfaces"] for m in (i["motif"] or {}).get("motifs", [])])]
    path = outdir / "scrims.cxc"
    path.write_text("\n".join(L) + "\n")
    return path
