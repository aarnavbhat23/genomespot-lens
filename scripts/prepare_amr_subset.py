"""Normalize a selected ATB subset from the observed AMRFinderPlus v4.2.5 schema.

The primary ATB result header was checked directly on 20 October 2026:
Name, Protein id, Contig id, Start, Stop, Strand, Element symbol, ... Type, Class.
Results and status files must come from the same completed release.
"""
import argparse
import json
from pathlib import Path

import pandas as pd

from genomespot_lens.amp_study import require_columns, assert_unique
from genomespot_lens.core import sha256


def normalize_chunk(chunk, ids, core_only=False):
    columns = ["Name", "Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol", "Type", "Class"]
    require_columns(chunk, columns, "AMRFinderPlus v4.2.5 results")
    selected = chunk[chunk.Name.isin(ids)].copy()
    if core_only:
        require_columns(selected, ["Scope"], "AMR core-only filter")
        selected = selected[selected.Scope.str.lower().eq("core")]
    locus = selected[["Contig id", "Start", "Stop", "Strand"]].astype(str).agg(":".join, axis=1)
    protein = selected["Protein id"].fillna("-").astype(str)
    locus = protein.where(~protein.isin(["-", "", "NA"]), locus)
    return pd.DataFrame({"genome_id": selected.Name, "element_type": selected.Type,
                         "gene_symbol": selected["Element symbol"], "drug_class": selected.Class,
                         "locus_id": locus})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True)
    parser.add_argument("--status", required=True)
    parser.add_argument("--genomes", required=True, help="Canonical genomes.csv; genome_id must be the ATB sample accession")
    parser.add_argument("--out", required=True)
    parser.add_argument("--core-only", action="store_true")
    args = parser.parse_args()
    genomes = pd.read_csv(args.genomes)
    require_columns(genomes, ["genome_id"], "genomes")
    ids = set(genomes.genome_id)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    collected = []
    for chunk in pd.read_csv(args.results, sep="\t", chunksize=50000, low_memory=False):
        collected.append(normalize_chunk(chunk, ids, args.core_only))
    normalized = pd.concat(collected, ignore_index=True)
    normalized.to_csv(out / "amr_hits.csv", index=False)
    status_rows = []
    for chunk in pd.read_csv(args.status, sep="\t", chunksize=100000):
        require_columns(chunk, ["sample", "status"], "AMRFinderPlus status")
        status_rows.append(chunk[chunk["sample"].isin(ids)][["sample", "status"]])
    status = pd.concat(status_rows, ignore_index=True).rename(columns={"sample": "genome_id"})
    assert_unique(status, ["genome_id"], "AMR status")
    status.to_csv(out / "amr_status.csv", index=False)
    report = {"results_sha256": sha256(args.results), "status_sha256": sha256(args.status), "core_only": args.core_only,
              "requested_genomes": len(ids), "status_rows_found": len(status), "hit_rows_found": len(normalized),
              "missing_status_ids": sorted(ids - set(status.genome_id)), "zero_rule": "Only PASS plus absent hits becomes zero during cohort preparation."}
    (out / "amr_import.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "missing_status_ids"}, indent=2))


if __name__ == "__main__":
    main()
