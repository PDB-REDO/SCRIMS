"""
Identify the common bait chain across a set of AlphaBridge predictions.

Detection is based on exact amino-acid sequence identity — not chain letter,
not sequence length. A chain is the bait if and only if its full sequence
string is present in every prediction in the screen.
"""

import json
from pathlib import Path


def load_chain_sequences(ab_json: Path) -> dict[str, str]:
    """Return {label_asym_id: sequence} for all polymer chains in one alphabridge_data.json."""
    data = json.loads(ab_json.read_text())
    return {
        c["label_asym_id"]: c["sequence"]
        for c in data["structure"][0]["chains"]["polymer"]
    }


def detect_bait(screen_dir: str | Path) -> dict:
    """
    Scan all alphabridge_data.json files under screen_dir and return a dict
    describing the common bait chain:

        {
            "sequence":     "MDFT...",   # full amino-acid sequence (the stable identifier)
            "length":       1337,
            "n_predictions": 27,
            "chain_ids": {               # chain label per prediction folder
                "fold_2026_02_04_pds5a_wapl":     "A",
                "fold_2026_02_04_pds5a_sororin_s1": "A",
                ...
            }
        }

    Bait identification is purely by sequence identity: the sequence that
    appears in every prediction is the bait. Length is never used as a
    criterion.

    Raises ValueError if zero or more than one sequence is shared across
    all predictions (ambiguous or no bait found).
    """
    screen_dir = Path(screen_dir)

    ab_jsons = sorted(screen_dir.glob("*/AlphaBridge/alphabridge_data.json"))
    if not ab_jsons:
        raise FileNotFoundError(f"No alphabridge_data.json found under {screen_dir}")

    per_prediction: list[dict[str, str]] = [
        load_chain_sequences(p) for p in ab_jsons
    ]

    # Intersect sequence sets across all predictions — only the bait survives
    sequence_sets = [set(chains.values()) for chains in per_prediction]
    common_sequences = sequence_sets[0].intersection(*sequence_sets[1:])

    if len(common_sequences) == 0:
        raise ValueError(
            f"No sequence is shared across all {len(ab_jsons)} predictions. "
            "Cannot auto-detect bait."
        )
    if len(common_sequences) > 1:
        raise ValueError(
            f"{len(common_sequences)} sequences are shared across all predictions — "
            "cannot unambiguously identify the bait."
        )

    bait_sequence = next(iter(common_sequences))

    # Build per-prediction mapping: folder name → chain label for the bait
    chain_ids = {}
    for ab_json, chains in zip(ab_jsons, per_prediction):
        pred_name = ab_json.parts[-3]   # .../fold_xxx/AlphaBridge/alphabridge_data.json
        chain_ids[pred_name] = next(
            cid for cid, seq in chains.items() if seq == bait_sequence
        )

    return {
        "sequence":      bait_sequence,
        "length":        len(bait_sequence),
        "n_predictions": len(ab_jsons),
        "chain_ids":     chain_ids,
    }


def detect_bait_from_partner_map(screen_dir: str | Path, partner_map: dict) -> dict:
    """
    Like detect_bait(), but sources bait-chain identity directly from a
    DB-recorded partner_map (prediction_id -> {"bait_chain": ..., ...}), as
    written by db_ingest_bait.py from the AF3 DB's own `chain` field -- the
    same field trusted over the CIF header's own (sometimes "poisoned")
    _struct_ref_seq block (see IMPROVEMENTS.md #15). Unlike detect_bait()'s
    sequence-intersection heuristic, this works even for a bait with only one
    prediction, since there's no other prediction needed to disambiguate against.

    Raises ValueError if the DB-recorded bait chain resolves to more than one
    distinct sequence across predictions -- a data inconsistency, not expected.
    """
    screen_dir = Path(screen_dir)
    ab_jsons = sorted(screen_dir.glob("*/AlphaBridge/alphabridge_data.json"))
    if not ab_jsons:
        raise FileNotFoundError(f"No alphabridge_data.json found under {screen_dir}")

    chain_ids = {}
    bait_sequences = set()
    for ab_json in ab_jsons:
        pred_name = ab_json.parts[-3]   # .../<prediction_id>/AlphaBridge/alphabridge_data.json
        bc = partner_map[pred_name]["bait_chain"]
        seq = load_chain_sequences(ab_json)[bc]
        chain_ids[pred_name] = bc
        bait_sequences.add(seq)

    if len(bait_sequences) > 1:
        raise ValueError(
            f"DB-recorded bait_chain resolves to {len(bait_sequences)} distinct sequences "
            f"across {len(ab_jsons)} predictions -- expected exactly 1."
        )

    bait_sequence = next(iter(bait_sequences))
    return {
        "sequence":      bait_sequence,
        "length":        len(bait_sequence),
        "n_predictions": len(ab_jsons),
        "chain_ids":     chain_ids,
    }


def detect_bait_from_structures(screen_dir: str | Path) -> dict:
    """`detect_bait`'s result without AlphaBridge: sequences come from each
    prediction's own structure file instead of alphabridge_data.json.

    Same rule, unchanged — a chain is the bait if and only if its full sequence
    appears in every prediction of the screen. Chain letter and length are never
    used as criteria. This is what lets AlphaSIFT run on raw AlphaFold output with
    no AlphaBridge run in between.
    """
    from alphasift.readers import discover, load_sequences

    screen = Path(screen_dir)
    per_pred = {}
    for name in discover(screen):
        seqs = load_sequences(screen, name)
        if seqs:
            per_pred[name] = seqs
    if not per_pred:
        raise FileNotFoundError(f"no readable predictions under {screen}")

    common = set.intersection(*(set(s.values()) for s in per_pred.values()))
    if not common:
        raise ValueError(
            f"no sequence is shared by all {len(per_pred)} predictions in {screen} — "
            "this does not look like a one-bait screen")
    if len(common) > 1:
        # Deterministic and explicit: the longest shared sequence, but say so.
        import warnings
        warnings.warn(f"{len(common)} sequences are common to every prediction in "
                      f"{screen.name}; taking the longest as the bait")
    bait_seq = max(common, key=len)
    chain_ids = {}
    for name, seqs in per_pred.items():
        chain_ids[name] = next(ch for ch, s in seqs.items() if s == bait_seq)
    return {"sequence": bait_seq, "length": len(bait_seq),
            "n_predictions": len(per_pred),
            "chain_ids": dict(sorted(chain_ids.items()))}


def bait_from_chain(screen_dir: str | Path, chain: str) -> dict:
    """The bait stated outright, as `{chain}` in every prediction.

    The sequence-intersection rule cannot resolve a bait from a SINGLE
    prediction — every chain is trivially common to all one of them — and a
    screen may legitimately contain two shared sequences. This is the explicit
    way out: name the chain.
    """
    from alphasift.readers import discover, load_sequences

    screen = Path(screen_dir)
    chain_ids, seq = {}, None
    for name in discover(screen):
        seqs = load_sequences(screen, name)
        if not seqs:
            continue
        if chain not in seqs:
            raise ValueError(f"prediction {name} has no chain {chain!r} "
                             f"(chains: {sorted(seqs)})")
        chain_ids[name] = chain
        if seq is None:
            seq = seqs[chain]
        elif seqs[chain] != seq:
            raise ValueError(f"chain {chain!r} of {name} differs in sequence from the "
                             "first prediction's — it cannot be a common bait")
    if not chain_ids:
        raise FileNotFoundError(f"no readable predictions under {screen}")
    return {"sequence": seq, "length": len(seq), "n_predictions": len(chain_ids),
            "chain_ids": dict(sorted(chain_ids.items()))}
