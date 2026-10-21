"""Reconstruct the archaeome input and source joins from checksum-pinned public files."""
import argparse
from collections import defaultdict, Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import re
from pathlib import Path

from genomespot_lens.archaeome_archive import (
    ARCHIVE_URL, ARCHIVE_MD5, CATALOG_SHA256, EXPECTED,
    read_archaeal_records, distinct_fragment_count, map_catalog,
)
from genomespot_lens.core import sha256


def write_csv(path, rows, fields=None):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--selected", type=Path, help="Optional published 265-row source table for independent attribution validation")
    parser.add_argument("--genus-table", type=Path, help="Optional transcribed published genus table for independent aggregate validation")
    args = parser.parse_args()
    with args.archive.open("rb") as handle:
        archive_md5 = hashlib.file_digest(handle, "md5").hexdigest()
    if archive_md5 != ARCHIVE_MD5 or sha256(args.catalog) != CATALOG_SHA256:
        raise ValueError("Input checksum mismatch; do not relabel a different release as this reconstruction")
    if args.out.exists():
        raise FileExistsError(args.out)
    records = read_archaeal_records(args.archive)
    with args.catalog.open() as handle:
        reader = csv.DictReader(handle)
        sequence_column = reader.fieldnames[0]
        catalog = list(reader)
    sequences = [r[sequence_column] for r in catalog]
    observed = {"protein_records": len(records),
                "distinct_proteins": len({r["sequence"] for r in records}),
                "source_taxa": len({r["taxid"] for r in records}),
                "distinct_screened_fragments": distinct_fragment_count(r["sequence"] for r in records),
                "positive_catalog_rows": len(sequences)}
    if observed != EXPECTED:
        raise ValueError(f"Published-count validation failed: {observed}")
    mappings = map_catalog(records, sequences)
    if len({r["catalog_row"] for r in mappings}) != len(sequences):
        raise ValueError("Some published positives do not map to the reconstructed source proteins")
    by_taxon, positive_by_taxon, taxa_by_peptide = defaultdict(list), defaultdict(set), defaultdict(set)
    for record in records:
        by_taxon[record["taxid"]].append(record)
    for row in mappings:
        positive_by_taxon[row["taxid"]].add(row["catalog_row"])
        taxa_by_peptide[row["catalog_row"]].add(row["taxid"])
        row["sequence_sha256"] = hashlib.sha256(sequences[row["catalog_row"]].encode()).hexdigest()
    sources = []
    for taxid, group in sorted(by_taxon.items()):
        sources.append({"taxon_id": taxid, "organism": group[0]["organism"],
            "reviewed_protein_records": len(group),
            "distinct_reviewed_proteins": len({r["sequence"] for r in group}),
            "distinct_screened_fragments": distinct_fragment_count(r["sequence"] for r in group),
            "distinct_high_scoring_peptides": len(positive_by_taxon[taxid])})
    selected_validation = None
    if args.selected:
        if sha256(args.selected) != "51ad7cc136604d56632bad87f3d2efae5adfcaa7513d922575c62d1b92953c7a":
            raise ValueError("Selected source table checksum mismatch")
        lookup = {r["accession"]: r for r in records}
        with args.selected.open() as handle:
            selected_reader = csv.DictReader(handle)
            selected_sequence_column = selected_reader.fieldnames[0]
            selected = list(selected_reader)
        mentions = 0
        for row in selected:
            taxa = set(re.findall(r"OX=(\d+)", row["ProteinInfo"]))
            for accession in re.findall(r"sp\|([^|]+)\|", row["ProteinInfo"]):
                mentions += 1
                if (accession not in lookup or row[selected_sequence_column] not in lookup[accession]["sequence"]
                        or lookup[accession]["taxid"] not in taxa):
                    raise ValueError("Published selected source attribution failed reconstruction validation")
        if len(selected) != 265 or mentions != 285:
            raise ValueError("Unexpected selected source attribution counts")
        selected_validation = {"table_sha256": sha256(args.selected), "rows": len(selected),
                               "validated_parent_accession_sequence_taxon_mentions": mentions}
    genus_validation = None
    if args.genus_table:
        counts, coverage = Counter(), Counter()
        aliases = {"Lokiarchaeum": "Candidatus Lokiarchaeum"}
        for row in sources:
            tokens = row["organism"].split()
            genus = " ".join(tokens[:2]) if tokens[0] == "Candidatus" else tokens[0]
            genus = aliases.get(genus, genus)
            counts[genus] += row["distinct_high_scoring_peptides"]
            coverage[genus] += 1
        with args.genus_table.open() as handle:
            published = list(csv.DictReader(handle))
        for row in published:
            genus = row["genus_as_published"]
            if coverage[genus] == 0 or counts[genus] != int(row["active_eps"]):
                raise ValueError(f"Published genus aggregate validation failed: {genus}")
        genus_validation = {"table_sha256": sha256(args.genus_table), "matched_genera": len(published),
            "aggregation": "Sum of distinct per-taxon positive counts, not genus-wide sequence deduplication",
            "genus_label_aliases": aliases,
            "caveat": "Published top/bottom genera are outcome-selected validation context, not an inference cohort"}
    args.out.mkdir(parents=True)
    write_csv(args.out / "source_taxa.csv", sources)
    write_csv(args.out / "positive_parent_mappings.csv", sorted(mappings,
              key=lambda r: (r["catalog_row"], r["taxid"], r["accession"], r["start"])),
              ["catalog_row", "sequence_sha256", "accession", "taxid", "start", "end"])
    write_csv(args.out / "protein_registry.csv", [{
        "accession": r["accession"], "taxon_id": r["taxid"], "isoform": r["isoform"],
        "length_aa": len(r["sequence"]), "sequence_sha256": hashlib.sha256(r["sequence"].encode()).hexdigest(),
        "refseq_cross_references": json.dumps(r["refseq"]),
        "proteome_cross_references": json.dumps(r["proteomes"])} for r in records])
    # Local FASTA is reproducible but not needed in the compact committed audit.
    with (args.out / "reviewed_proteins.fasta").open("w") as handle:
        for r in records:
            handle.write(f'>{r["accession"]} OX={r["taxid"]}\n{r["sequence"]}\n')
    audit = {"created_utc": datetime.now(timezone.utc).isoformat(),
        "description": "Count-validated UniProt 2023_03 reconstruction; not an author-supplied accession manifest",
        "paper": "https://doi.org/10.1038/s41564-025-02061-0",
        "paper_uniprot_access_date": "2023-08-24", "archive_url": ARCHIVE_URL,
        "archive_bytes": args.archive.stat().st_size, "archive_md5": archive_md5,
        "archive_sha256": sha256(args.archive), "catalog_sha256": sha256(args.catalog),
        "expected_counts": EXPECTED, "observed_counts": observed,
        "all_five_counts_match": observed == EXPECTED,
        "selected_source_validation": selected_validation,
        "published_genus_validation": genus_validation,
        "all_positive_peptides_mapped": True, "parent_occurrences": len(mappings),
        "multi_taxon_positive_peptides": sum(len(v) > 1 for v in taxa_by_peptide.values()),
        "taxa_with_positive_hits": sum(bool(positive_by_taxon[t]) for t in by_taxon),
        "taxa_with_zero_catalog_hits": sum(not positive_by_taxon[t] for t in by_taxon),
        "endpoint": "distinct published mean-panel MIC <=100 micromolar positive sequences per source taxon",
        "exposure": "distinct canonical 8-50 aa substrings in each taxon's reviewed source proteins, including isoforms",
        "unit": "source taxon, not a resolved whole-genome assembly",
        "caveat": "Whole-genome mapping, quality, taxonomy/tree and independent habitat remain required. Reviewed-protein coverage is highly uneven. Multiple-taxon mappings are retained, not forced to one source. Published predictions were not rerun and activity was not measured here.",
        "license": "UniProt CC BY 4.0; retain upstream paper and UniProt attribution",
        "outputs": {p.name: {"sha256": sha256(p), "bytes": p.stat().st_size}
                    for p in sorted(args.out.iterdir())}}
    (args.out / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({k: audit[k] for k in ["observed_counts", "all_positive_peptides_mapped", "parent_occurrences", "multi_taxon_positive_peptides", "taxa_with_positive_hits", "taxa_with_zero_catalog_hits"]}, indent=2))


if __name__ == "__main__":
    main()
