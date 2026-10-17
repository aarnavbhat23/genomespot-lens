from collections import Counter
import random

import pytest

from genomespot_lens.core import map_sequences, matched_sample, numeric_prediction, quality_flags, read_fasta, validate_ids
from genomespot_lens.cli import summarize


@pytest.mark.parametrize("bad", ["", "ACTG", ">a\n", ">a\nACD\n>a\nAAA", ">\nACD", ">a\nAC*D"])
def test_bad_fasta_rejected(tmp_path, bad):
    path = tmp_path / "bad.faa"
    path.write_text(bad)
    with pytest.raises(ValueError):
        read_fasta(path)


def test_fasta_terminal_stop_and_case(tmp_path):
    path = tmp_path / "good.faa"
    path.write_text(">a description\nac\nd*\n")
    assert read_fasta(path) == {"a": "ACD"}


def test_exact_mapping_does_not_guess_paralogs_or_partial_matches():
    ids, records = map_sequences({"a": "ACD", "b": "ACD", "c": "FGH"}, {"x": "ACD", "y": "FGH", "z": "FG"})
    assert ids == ["c"]
    assert [r["status"] for r in records] == ["ambiguous", "matched", "unmatched"]


def test_matched_controls_preserve_count_length_bins_and_reproducibility():
    lengths = {str(i): 30 + i * 11 for i in range(100)}
    subset = list(lengths)[::4]
    control = matched_sample(lengths, subset, random.Random(23), 100)
    assert control == matched_sample(lengths, subset, random.Random(23), 100)
    assert len(control) == len(set(control)) == len(subset)
    assert Counter(lengths[k] // 100 for k in control) == Counter(lengths[k] // 100 for k in subset)
    assert set(control) <= lengths.keys()
    assert set(control) != set(subset)


@pytest.mark.parametrize("label,prob,want", [("tolerant", .9, .9), ("not tolerant", .8, .2), ("not intolerant", .8, .2)])
def test_probability_is_probability_of_tolerance(label, prob, want):
    assert numeric_prediction("oxygen", {"value": label, "error": prob}) == pytest.approx(want)


def test_quality_separates_benign_salinity_clamp_and_novelty():
    p = {"value": 0, "warning": "min_exceeded", "is_novel": False}
    assert quality_flags("salinity_min", p) == []
    assert quality_flags("temperature_optimum", p) == ["min_exceeded"]
    p["is_novel"] = True
    assert quality_flags("salinity_min", p) == ["outside_training_distribution"]


def test_unknown_and_duplicate_ids_rejected(tmp_path):
    path = tmp_path / "ids.txt"
    for text in ["z\n", "a\na\n", ""]:
        path.write_text(text)
        with pytest.raises(ValueError):
            validate_ids(path, {"a": "ACD"})


def test_missing_predictions_are_counted_not_converted_to_zero():
    def p(value, novel=False):
        return {"temperature_optimum": {"value": value, "is_novel": novel, "warning": None}}
    records = [{"kind": "full", "predictions": p(30)}, {"kind": "selected", "predictions": p(35)}, {"kind": "count_matched", "predictions": p(None)}, {"kind": "count_matched", "predictions": p(32, True)}]
    summary = summarize(records, ["temperature_optimum"])[0]
    assert summary["control_total"] == 2
    assert summary["control_finite"] == 1
    assert summary["control_unflagged"] == 0
    assert summary["control_delta_median"] == 2
    assert summary["delta"] == 5


def test_oxygen_low_confidence_and_missing():
    assert "low_class_confidence" in quality_flags("oxygen", {"value": "tolerant", "error": .75, "is_novel": False})
    assert numeric_prediction("oxygen", {"value": "tolerant", "error": 0.2}) is None
