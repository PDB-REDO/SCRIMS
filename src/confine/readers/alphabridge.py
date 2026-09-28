"""
AlphaBridge output as a second confidence source (AlphaFold3 only).

AlphaFold2 predictions expose PAE and nothing else usable here, so the PAE
channel is the only option for them. AlphaFold3 additionally carries a
distogram, from which AlphaBridge derives contact probabilities and sweeps a
cut-off from 0 to 1, recording at each level which residue ranges form an
interface. The value taken per bait residue is the **highest cut-off at which
that residue still belongs to an interface** — a model-derived confidence in the
contact, not a distance proxy, which is why it is offered alongside the PAE
channel rather than as a legacy fallback.

The two sources share the whole downstream pipeline but not their scales: on the
PAE channel 0.8 means PAE 5 A, whereas here it means the residue survives to
cut-off 0.8 in AlphaBridge's sweep. Threshold parameters keep their numeric
defaults across sources but not their interpretation.

`bait_row` is extracted from viz_int.patches.load_screen_matrix.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def ab_json(screen_dir, name: str) -> Path | None:
    """The AlphaBridge JSON for one prediction, if the run produced one."""
    p = Path(screen_dir) / name / "AlphaBridge" / "alphabridge_data.json"
    return p if p.exists() else None


def bait_row(path, bait_chain: str, length: int) -> np.ndarray:
    """Per bait residue, the maximum sweep cut-off at which it is in an interface.

    Zero means the residue never appears in an interface at any cut-off, which
    is the same "absent" encoding the PAE channel uses, so the statistics run
    unchanged.
    """
    d = json.loads(Path(path).read_text())
    v = np.zeros(length)
    for e in d["interactions"]:
        c = e["cut-off"]
        for iface in e["interfaces"]:
            for lk in iface["links"]:
                for end in ("first", "second"):
                    if lk[end]["asym_id"] == bait_chain:
                        s, en = lk[end]["link_range"]["start"], lk[end]["link_range"]["end"]
                        a, b = max(0, s - 1), min(length, en)
                        v[a:b] = np.maximum(v[a:b], c)
    return v
