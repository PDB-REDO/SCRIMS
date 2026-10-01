"""A synthetic AF3-style screen, so the tests need no real prediction data."""
import json
import numpy as np
import pytest

COLS = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_comp_id",
        "label_asym_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z",
        "B_iso_or_equiv"]


def _cif(bait_len, partner_len, contact_residues, partner_aa="LEU", jitter=0.0):
    """Bait along x at y=0; partner along x at y=20 except where it contacts the
    bait, where it drops to y=4 (inside the 8 A Cb selection bar)."""
    rows, aid = [], 1
    for i in range(1, bait_len + 1):
        for atom, dz in (("N", -0.5), ("CA", 0.0), ("C", 0.5), ("CB", 1.2)):
            rows.append(["ATOM", aid, atom[0], atom, "ALA", "A", i,
                         3.5 * i, 0.0 + jitter, dz, 90.0]); aid += 1
    for j in range(1, partner_len + 1):
        y = 4.0 if j in contact_residues else 20.0
        for atom, dz in (("N", -0.5), ("CA", 0.0), ("C", 0.5), ("CB", 1.2)):
            rows.append(["ATOM", aid, atom[0], atom, partner_aa, "B", j,
                         3.5 * contact_residues.get(j, j), y, dz, 85.0]); aid += 1
    out = ["data_test", "#", "loop_"] + [f"_atom_site.{c}" for c in COLS]
    out += [" ".join(str(v) for v in r) for r in rows] + ["#"]
    return "\n".join(out)


def _conf(bait_len, partner_len, good_pae, contact_residues):
    T = bait_len + partner_len
    pae = np.full((T, T), 28.0)
    np.fill_diagonal(pae, 0.3)
    for j, bait_res in contact_residues.items():
        pj = bait_len + j - 1
        bi = bait_res - 1
        pae[bi, pj] = pae[pj, bi] = good_pae
    return {"pae": pae.tolist(),
            "token_chain_ids": ["A"] * bait_len + ["B"] * partner_len,
            "token_res_ids": list(range(1, bait_len + 1)) + list(range(1, partner_len + 1))}


@pytest.fixture
def screen(tmp_path):
    """One 40-residue bait, 6 partners. Four of them confidently engage residues
    10-14 (the planted convergent site); two engage scattered positions."""
    bait_len = 40
    plan = {
        "p1": ({1: 10, 2: 11, 3: 12, 4: 13, 5: 14}, 1.0),
        "p2": ({1: 10, 2: 11, 3: 12, 4: 13, 5: 14}, 1.5),
        "p3": ({1: 11, 2: 12, 3: 13, 4: 14, 5: 15}, 2.0),
        "p4": ({1: 10, 2: 11, 3: 12, 4: 13, 5: 14}, 2.5),
        "p5": ({1: 30, 2: 31, 3: 32, 4: 33, 5: 34}, 9.5),
        "p6": ({1: 2, 2: 20, 3: 36}, 9.8),
    }
    # Each partner gets its own residue type, so the bait is the ONLY sequence
    # shared by every prediction -- which is what bait detection must key on.
    aa = {"p1": "LEU", "p2": "ILE", "p3": "VAL", "p4": "PHE", "p5": "TRP", "p6": "TYR"}
    for name, (contacts, pae) in plan.items():
        d = tmp_path / name
        d.mkdir()
        (d / f"{name}_model.cif").write_text(_cif(bait_len, 6, contacts, aa[name]))
        (d / f"{name}_confidences.json").write_text(
            json.dumps(_conf(bait_len, 6, pae, contacts)))
    return tmp_path


AA1 = {"p1": "L", "p2": "I", "p3": "V", "p4": "F", "p5": "W", "p6": "Y"}


def add_alphabridge(screen, plan=None):
    """Attach a synthetic AlphaBridge run to predictions of a screen.

    `plan` maps prediction name -> (bait_start, bait_end, highest cut-off at
    which that range is still in an interface). Written in the shape
    alphasift.readers.alphabridge.bait_row reads: a sweep of cut-offs, each
    listing the interfaces surviving at that level.
    """
    plan = plan or {"p1": (10, 14, 0.9), "p2": (10, 14, 0.7), "p3": (11, 15, 0.5)}
    for name, (a, b, top) in plan.items():
        d = screen / name / "AlphaBridge"
        d.mkdir(parents=True, exist_ok=True)
        interactions = []
        for c in [round(0.1 * k, 1) for k in range(1, 11)]:
            if c > top:
                break
            interactions.append({
                "cut-off": c,
                "interfaces": [{"links": [{
                    "first": {"asym_id": "A", "link_range": {"start": a, "end": b}},
                    "second": {"asym_id": "B", "link_range": {"start": 1, "end": 5}},
                }]}],
            })
        # Real AlphaBridge output carries the chain sequences too; bait detection
        # prefers them when present.
        doc = {"structure": [{"chains": {"polymer": [
                   {"label_asym_id": "A", "sequence": "A" * 40},
                   {"label_asym_id": "B", "sequence": AA1[name] * 6}]}}],
               "interactions": interactions}
        (d / "alphabridge_data.json").write_text(json.dumps(doc))
    return screen
