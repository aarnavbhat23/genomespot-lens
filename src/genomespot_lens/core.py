"""FASTA validation, subset mapping, matched controls and interpretation."""

from collections import Counter, defaultdict
import gzip
import hashlib
import math
from pathlib import Path
import random


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_fasta(path, protein=True):
    """Read unique first-token IDs; reject empty, malformed or duplicate records."""
    opener = gzip.open if str(path).endswith(".gz") else open
    records, name, parts = {}, None, []
    alphabet = set("ACDEFGHIKLMNPQRSTVWYBXZJUO" if protein else "ACGTRYSWKMBDHVN")

    def save():
        if name is None:
            return
        seq = "".join(parts).upper()
        if protein:
            seq = seq.removesuffix("*")
        if not seq or set(seq) - alphabet:
            raise ValueError(f"Empty or invalid sequence: {name}")
        if name in records:
            raise ValueError(f"Duplicate FASTA ID: {name}")
        records[name] = seq

    with opener(path, "rt") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                if not line[1:].split():
                    raise ValueError("Empty FASTA header")
                name, parts = line[1:].split()[0], []
            else:
                if name is None:
                    raise ValueError("Sequence before FASTA header")
                parts.append(line)
    save()
    if not records:
        raise ValueError("Empty FASTA")
    return records


def map_sequences(full, subset):
    """Exact full-length matching only. Ambiguous gene copies remain unresolved."""
    by_sequence = defaultdict(list)
    for identifier, sequence in full.items():
        by_sequence[sequence].append(identifier)
    matched, details = set(), []
    for identifier, sequence in subset.items():
        candidates = by_sequence.get(sequence, [])
        status = "matched" if len(candidates) == 1 else "ambiguous" if candidates else "unmatched"
        if status == "matched":
            matched.add(candidates[0])
        details.append({"source_id": identifier, "status": status, "candidates": candidates})
    return sorted(matched), details


def matched_sample(lengths, selected, rng, bin_width=None):
    """Draw from the full proteome, including selected IDs, without replacement.

    Width W matches counts in fixed bins [0,W), [W,2W), ...; not exact lengths.
    """
    if not selected or not set(selected) <= lengths.keys():
        raise ValueError("Subset must be nonempty and contained in the full proteome")
    if bin_width is None:
        return sorted(rng.sample(sorted(lengths), len(selected)))
    if bin_width < 1:
        raise ValueError("Length bin width must be positive")
    pools = defaultdict(list)
    for identifier in sorted(lengths):
        pools[lengths[identifier] // bin_width].append(identifier)
    counts = Counter(lengths[x] // bin_width for x in selected)
    return sorted(x for bucket, count in sorted(counts.items()) for x in rng.sample(pools[bucket], count))


def numeric_prediction(target, prediction):
    value, error = prediction.get("value"), prediction.get("error")
    if value is None:
        return None
    if target == "oxygen":
        if error is None or not math.isfinite(float(error)) or not 0.5 <= float(error) <= 1:
            return None
        if value == "tolerant":
            return float(error)
        if value in {"not tolerant", "not intolerant"}:
            return 1 - float(error)
        raise ValueError(f"Unknown oxygen label: {value}")
    value = float(value)
    return value if math.isfinite(value) else None


def quality_flags(target, prediction):
    flags = []
    if numeric_prediction(target, prediction) is None:
        flags.append("missing_prediction")
    if prediction.get("is_novel") is True:
        flags.append("outside_training_distribution")
    elif prediction.get("is_novel") is None:
        flags.append("novelty_unknown")
    warning = prediction.get("warning")
    benign_zero = target in {"salinity_min", "salinity_optimum"} and prediction.get("value") == 0 and warning == "min_exceeded"
    if warning and not benign_zero:
        flags.append(str(warning))
    if target == "oxygen" and prediction.get("error") is not None and float(prediction["error"]) <= 0.75:
        flags.append("low_class_confidence")
    return flags


def validate_ids(path, records):
    ids = [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("Subset IDs must be nonempty and unique")
    unknown = set(ids) - records.keys()
    if unknown:
        raise ValueError(f"Unknown subset IDs: {sorted(unknown)[:5]}")
    return sorted(ids)
