"""Exact fragment accounting must agree with brute force, especially at boundaries."""
import random
import pytest

pytest.importorskip("ahocorasick")
pytest.importorskip("pydivsufsort")
from genomespot_lens.archaeome_archive import distinct_fragment_count, map_catalog


def brute(sequences, lo, hi):
    canonical = set("ACDEFGHIKLMNPQRSTVWY")
    return len({s[i:i+k] for s in sequences for k in range(lo, hi+1)
                for i in range(len(s)-k+1) if set(s[i:i+k]) <= canonical})


def test_fragment_count_matches_exhaustive_sets():
    rng = random.Random(23)
    for _ in range(40):
        sequences = ["".join(rng.choices("AACDXY", k=rng.randrange(0, 35))) for _ in range(6)]
        assert distinct_fragment_count(sequences, 2, 8) == brute(sequences, 2, 8)
    assert distinct_fragment_count(["ACDEFGHIK", "ACDEFGHIK", "ACDEFXGHIK"], 8, 50) == 3
    assert distinct_fragment_count([], 8, 50) == 0
    with pytest.raises(ValueError):
        distinct_fragment_count(["ACD"], 5, 3)


def test_mapping_retains_multiple_taxa_and_overlapping_occurrences():
    records = [{"accession": "a", "taxid": "1", "sequence": "AAAA"},
               {"accession": "b", "taxid": "2", "sequence": "CAAA"}]
    mapped = map_catalog(records, ["AAA", "CA"])
    assert len(mapped) == 4
    assert [(r["start"], r["end"]) for r in mapped if r["accession"] == "a"] == [(1, 3), (2, 4)]
    assert {r["taxid"] for r in mapped if r["catalog_row"] == 0} == {"1", "2"}
    with pytest.raises(ValueError):
        map_catalog(records, ["AAA", "AAA"])


def test_committed_reconstruction_has_complete_consistent_source_links():
    import csv
    import hashlib
    import json
    from collections import defaultdict
    from pathlib import Path
    root = Path(__file__).resolve().parents[1] / "examples/amp-study/archaeome-reconstruction"
    audit = json.loads((root / "audit.json").read_text())
    for name in ["protein_registry.csv", "positive_parent_mappings.csv", "source_taxa.csv"]:
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == audit["outputs"][name]["sha256"]
    proteins = list(csv.DictReader((root / "protein_registry.csv").open()))
    assert len(proteins) == 19710
    assert len({r["sequence_sha256"] for r in proteins}) == 18677
    source = list(csv.DictReader((root / "source_taxa.csv").open()))
    assert len(source) == 233
    lookup = {r["accession"]: r for r in proteins}
    mapped = list(csv.DictReader((root / "positive_parent_mappings.csv").open()))
    assert {int(r["catalog_row"]) for r in mapped} == set(range(12623))
    by_taxon = defaultdict(set)
    for row in mapped:
        assert row["taxid"] == lookup[row["accession"]]["taxon_id"]
        assert 1 <= int(row["start"]) <= int(row["end"]) <= int(lookup[row["accession"]]["length_aa"])
        by_taxon[row["taxid"]].add(row["sequence_sha256"])
    for row in source:
        assert len(by_taxon[row["taxon_id"]]) == int(row["distinct_high_scoring_peptides"])
        assert int(row["distinct_screened_fragments"]) > 0


def test_real_batch_outputs_match_input_identity_and_prediction_provenance():
    import csv
    import hashlib
    import json
    from pathlib import Path
    root = Path(__file__).resolve().parents[1] / "examples/amp-study/reconstructed-genomes"
    audit = json.loads((root / "prediction_audit.json").read_text())
    identities = {r["assembly_accession"]: r for r in csv.DictReader((root / "source_protein_identity.csv").open())}
    rows = list(csv.DictReader((root / "genomespot.csv").open()))
    assert len(rows) == audit["prediction_rows"] == 1380
    assert len({(r["genome_id"], r["target"]) for r in rows}) == len(rows)
    assert set(identities) == {r["genome_id"] for r in rows}
    for name, expected_hash in audit["prediction_records"].items():
        path = root / "prediction_records" / name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_hash
        record = json.loads(path.read_text())
        identity = identities[record["genome_id"]]
        assert record["signature"]["proteins_sha256"] == identity["proteins_sha256"]
        assert record["signature"]["contigs_sha256"] == identity["contigs_sha256"]
        assert len(record["predictions"]) == 10
