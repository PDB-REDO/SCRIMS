"""
Input dispatch: one AlphaFold prediction, whatever produced it, reduced to the
four arrays the matrix needs.

Every reader returns a `Prediction` with per-token chain ids, residue ids, Cbeta
coordinates and a (T, T) PAE matrix, plus the per-chain sequences SCRIMS uses
to identify the bait. Adding a format means adding a reader here; nothing
downstream changes.

Formats:
  af3         AlphaFold3 server/local output: <pred>_model.cif + <pred>_confidences.json
  afdb        AlphaFold2 in AFDB/ColabFold bulk packaging:
              <id>-model_v1.cif[.zst] + <id>-predicted_aligned_error_v1.json[.zst]
  alphabridge an AlphaBridge run's alphabridge_data.json carrying a precomputed
              per-bait-residue reduction (for screens whose raw PAE is not kept
              on disk); handled in scrims.screen, since it bypasses the arrays

AlphaBridge is never required: it is one optional way in, not a dependency.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from scrims.readers.mmcif import parse_cif_cb, load_confidences
from scrims.readers.afdb import load_entry

AA3_TO_1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V',
}


@dataclass
class Prediction:
    """One prediction, format-independent."""
    name: str
    format: str
    chain_ids: np.ndarray     # (T,) per-token chain label
    res_ids: np.ndarray       # (T,) per-token residue number, per chain
    cb: np.ndarray            # (T, 3) Cbeta (pseudo-Cbeta for glycine)
    pae: np.ndarray           # (T, T) predicted aligned error, Angstrom
    res_names: np.ndarray | None = None   # (T,) 3-letter, for sequences

    def chain_sequences(self) -> dict[str, str]:
        """{chain label: 1-letter sequence}, in token order. Used to identify the
        bait by exact sequence identity when no AlphaBridge output is present."""
        if self.res_names is None:
            return {}
        seqs: dict[str, list[str]] = {}
        for ch, rn in zip(self.chain_ids, self.res_names):
            seqs.setdefault(str(ch), []).append(AA3_TO_1.get(str(rn).upper(), 'X'))
        return {k: ''.join(v) for k, v in seqs.items()}


def _read_af3(screen: Path, name: str) -> Prediction | None:
    d = screen / name
    cifs = sorted(d.glob("*_model.cif"))
    confs = [p for p in d.glob("*_confidences.json") if "summary" not in p.name]
    if not (cifs and confs):
        return None
    res = parse_cif_cb(cifs[0])
    conf = load_confidences(confs[0])
    return Prediction(name=name, format="af3", chain_ids=res.chain_ids,
                      res_ids=res.res_ids, cb=res.cb, pae=conf.pae,
                      res_names=res.res_names)


def _read_afdb(screen: Path, name: str) -> Prediction | None:
    for cif_suffix, pae_suffix in (("-model_v1.cif.zst", "-predicted_aligned_error_v1.json.zst"),
                                   ("-model_v1.cif", "-predicted_aligned_error_v1.json")):
        cif, pae = screen / f"{name}{cif_suffix}", screen / f"{name}{pae_suffix}"
        if cif.exists() and pae.exists():
            e = load_entry(str(cif), str(pae))
            return Prediction(name=name, format="afdb", chain_ids=e.token_chain_ids,
                              res_ids=e.token_res_ids, cb=e.cb, pae=e.pae,
                              res_names=e.comp_ids)
    return None


READERS = (_read_af3, _read_afdb)


def load_prediction(screen_dir, name: str) -> Prediction | None:
    """First reader that recognises `name` in `screen_dir`, else None."""
    screen = Path(screen_dir)
    for reader in READERS:
        got = reader(screen, name)
        if got is not None:
            return got
    return None


def discover(screen_dir) -> list[str]:
    """Prediction names in a screen directory, from whichever layout is present.

    A subdirectory counts if a reader recognises it; a loose
    `<name>-model_v1.cif[.zst]` counts as `<name>`. Names are returned sorted and
    de-duplicated, so a screen that has both (AFDB dumps keep the compressed
    files beside an unpacked folder) yields each prediction once.
    """
    screen = Path(screen_dir)
    names = set()
    for p in screen.iterdir():
        if p.is_dir() and not p.name.startswith(('_', '.')):
            names.add(p.name)
        elif "-model_v1.cif" in p.name:
            names.add(p.name.split("-model_v1.cif")[0])
    return sorted(names)


def load_sequences(screen_dir, name: str) -> dict[str, str] | None:
    """{chain label: sequence} for one prediction, parsing the structure file only.

    Bait detection needs sequences from every prediction in a screen but no PAE;
    the PAE matrices are the expensive part (and, for a large screen, far too
    large to hold all at once), so this exists as a cheap first pass.
    """
    screen = Path(screen_dir)
    d = screen / name
    cifs = sorted(d.glob("*_model.cif"))
    if cifs:
        res = parse_cif_cb(cifs[0])
        return Prediction(name, "af3", res.chain_ids, res.res_ids, res.cb,
                          np.zeros((0, 0)), res.res_names).chain_sequences()
    for suffix in ("-model_v1.cif.zst", "-model_v1.cif"):
        cif = screen / f"{name}{suffix}"
        if cif.exists():
            from scrims.readers.afdb import parse_cif, read_maybe_zst
            chain_ids, res_ids, comps, cb, _plddt = parse_cif(read_maybe_zst(str(cif)).decode())
            return Prediction(name, "afdb", chain_ids, res_ids, cb,
                              np.zeros((0, 0)), comps).chain_sequences()
    return None


def model_cif(screen_dir, name: str, tmpdir) -> Path:
    """A plain-text model CIF for `name`, for tools that need a file on disk (Bio.PDB,
    Foldseek). An existing `<name>/*_model.cif` or loose `<name>-model_v1.cif`
    is returned as is; a loose `.cif.zst` is decompressed into `tmpdir`."""
    screen = Path(screen_dir)
    cifs = sorted((screen / name).glob("*_model.cif"))
    if cifs:
        return cifs[0]
    if (screen / f"{name}-model_v1.cif").exists():
        return screen / f"{name}-model_v1.cif"
    zst = screen / f"{name}-model_v1.cif.zst"
    if zst.exists():
        from scrims.readers.afdb import read_maybe_zst
        out = Path(tmpdir) / f"{name}_model.cif"
        out.write_bytes(read_maybe_zst(str(zst)))
        return out
    raise FileNotFoundError(f"no structure file for {name} in {screen}")


def load_atoms(screen_dir, name: str, chain: str) -> dict:
    """{residue number: (n, 3) heavy atoms} for one chain of one prediction.

    Needed only by the 3-D interface merge, which works on the bait chain.
    """
    from scrims.readers.mmcif import parse_cif_atoms

    screen = Path(screen_dir)
    cifs = sorted((screen / name).glob("*_model.cif"))
    if cifs:
        return parse_cif_atoms(cifs[0], chain=chain)
    for suffix in ("-model_v1.cif.zst", "-model_v1.cif"):
        cif = screen / f"{name}{suffix}"
        if cif.exists():
            from scrims.readers.afdb import read_maybe_zst
            return parse_cif_atoms(cif, chain=chain,
                                   text=read_maybe_zst(str(cif)).decode())
    raise FileNotFoundError(f"no structure file for {name} in {screen}")
