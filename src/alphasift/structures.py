"""
Partner-side structure loading for the partner annotation and the motif engine.

For every binder prediction: the heavy atoms of the bait and partner chains, the
partner's residue names (PTMs kept, as their parent residue in lowercase) and,
for the motif engine's touch rule, one Cbeta per residue from the same parse.

`T2O`, `aa1`, `load_binder_structures` and `_cb_of_chain` are extracted verbatim
from viz_int.structures (2026-10-02), with one change: the model file of each
prediction is passed in (`cif_of`), because AlphaSIFT also reads screens whose
models exist only as loose `.cif.zst` (see alphasift.readers.model_cif). The
parse is Bio.PDB's, as in the source, so residue numbering and the
standard-residue filter are the validated ones.
"""

from __future__ import annotations

import numpy as np
from Bio.PDB import MMCIFParser

from alphasift.readers.mmcif import pseudo_cb

# three-letter -> one-letter, with the PTMs the screen actually contains
# (mapped to their parent residue in lowercase so they stay visible).
T2O = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEP': 's', 'TPO': 't', 'PTR': 'y', 'ALY': 'k', 'MLY': 'k', 'M3L': 'k', 'HYP': 'p', 'CSO': 'c',
}


def aa1(rn):
    """Three-letter residue name -> one-letter code ('X' if unknown)."""
    return T2O.get(rn, 'X')


def load_binder_structures(cif_of, bait, binder_preds, with_cb=False):
    """Load heavy-atom coordinates for each binder prediction.

    `cif_of` = {prediction: plain-text model CIF path}. Returns three dicts keyed by
    prediction name:
      bait_atoms[pr]  : {bait_residue_number   -> (n_atoms, 3) heavy-atom coords}
      int_atoms[pr]   : {interactor_residue_number -> (n_atoms, 3) heavy-atom coords}
      int_seq[pr]     : {interactor_residue_number -> three-letter residue name}

    The bait chain is `bait['chain_ids'][pr]`; the interactor chain is the other
    polymer chain in the model. Only standard (hetflag ' ') residues are kept.

    with_cb=True also returns, from the same parse,
      bait_cb[pr] / int_cb[pr] : {residue_number -> Cb xyz (float64)}
    with the virtual Cb from N/CA/C where Cb is missing (Gly) -- the construction of the
    confidence matrix; residues without CA are left out.
    """
    parser = MMCIFParser(QUIET=True)
    bait_atoms, int_atoms, int_seq = {}, {}, {}
    bait_cb, int_cb = {}, {}
    for pr in binder_preds:
        bc = bait['chain_ids'][pr]
        model = parser.get_structure(pr, str(cif_of[pr]))[0]
        icid = [c.id for c in model if c.id != bc][0]
        bait_atoms[pr] = {
            r.id[1]: np.array([a.coord for a in r if not a.name.startswith('H')], float)
            for r in model[bc] if r.id[0] == ' '
        }
        int_atoms[pr] = {
            r.id[1]: np.array([a.coord for a in r if not a.name.startswith('H')], float)
            for r in model[icid] if r.id[0] == ' '
        }
        int_seq[pr] = {r.id[1]: r.resname for r in model[icid] if r.id[0] == ' '}
        if with_cb:
            bait_cb[pr] = _cb_of_chain(model[bc])
            int_cb[pr] = _cb_of_chain(model[icid])
    if with_cb:
        return bait_atoms, int_atoms, int_seq, bait_cb, int_cb
    return bait_atoms, int_atoms, int_seq


def _cb_of_chain(chain):
    """{residue_number -> Cb xyz} (virtual Cb from N/CA/C when Cb is missing; residues without CA
    skipped) for the standard residues of one chain."""
    out = {}
    for r in chain:
        if r.id[0] != ' ':
            continue
        at = {a.get_name(): a.coord.astype(float) for a in r}
        if 'CA' not in at:
            continue
        if 'CB' in at:
            out[r.id[1]] = at['CB']
        elif all(k in at for k in ('N', 'CA', 'C')):
            out[r.id[1]] = pseudo_cb(at['N'], at['CA'], at['C'])
    return out
