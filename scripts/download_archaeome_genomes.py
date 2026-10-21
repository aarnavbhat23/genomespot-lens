"""Download resolved current source assemblies and verify historical protein coverage.

No AMP scores are used to choose inputs. Technical GenomeSPOT runs may include
partial historical-protein matches; those mismatches must be excluded or audited
before joining historical richness to current genome covariates.
"""
import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

import requests

from genomespot_lens.core import read_fasta, sha256
from source_genome_pilot import download


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader();writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--reconstruction", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    inputs = args.out / "inputs"; inputs.mkdir()
    source_proteins = read_fasta(args.reconstruction / "reviewed_proteins.fasta")
    by_taxon = defaultdict(set)
    for row in csv.DictReader((args.reconstruction / "protein_registry.csv").open()):
        by_taxon[row["taxon_id"]].add(source_proteins[row["accession"]])
    mappings = list(csv.DictReader(args.mapping.open()))
    status, manifest, identity = [], [], []
    for row in mappings:
        if row["status"] != "resolved_current_exact_taxon":
            continue
        taxid, accession = row["taxon_id"], row["assembly_accession"]
        item = {"taxon_id": taxid, "assembly_accession": accession, "status": "FAIL", "reason": ""}
        try:
            fna, faa = inputs / f"{accession}.fna.gz", inputs / f"{accession}.faa.gz"
            download(row["contigs_url"], fna);download(row["proteins_url"], faa)
            proteins = read_fasta(faa);read_fasta(fna, protein=False)
            whole = set(proteins.values())
            exact = sum(s in whole for s in by_taxon[taxid])
            # Containment handles a historical partial protein within a longer current annotation.
            contained = sum(s in whole or any(s in p for p in whole) for s in by_taxon[taxid])
            identity.append({"taxon_id": taxid, "assembly_accession": accession,
                "historical_distinct_source_proteins": len(by_taxon[taxid]),
                "exact_whole_protein_matches": exact, "contained_source_protein_matches": contained,
                "all_historical_sources_contained": contained == len(by_taxon[taxid]),
                "current_protein_count": len(proteins), "current_proteome_aa": sum(map(len, proteins.values())),
                "contigs_sha256": sha256(fna), "proteins_sha256": sha256(faa)})
            manifest.append({"genome_id": accession, "contigs": str(fna.relative_to(args.out)),
                             "proteins": str(faa.relative_to(args.out))})
            item.update(status="PASS", reason="Input valid; historical protein coverage recorded separately")
        except (requests.RequestException, ValueError, KeyError, OSError, EOFError) as exc:
            item["reason"] = str(exc)
        status.append(item)
        write_csv(args.out / "download_status.csv", status, ["taxon_id", "assembly_accession", "status", "reason"])
        write_csv(args.out / "genome_inputs.csv", manifest, ["genome_id", "contigs", "proteins"])
        if identity:
            write_csv(args.out / "source_protein_identity.csv", identity, list(identity[0]))
        print(len(status), accession, item["status"], flush=True)
    (args.out / "audit.json").write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(), "mapping_sha256": sha256(args.mapping),
        "historical_registry_sha256": sha256(args.reconstruction / "protein_registry.csv"),
        "attempted_resolved_genomes": len(status), "download_pass": len(manifest),
        "all_historical_sources_contained": sum(x["all_historical_sources_contained"] for x in identity),
        "scope": "Current exact-taxon assembly inputs; technical predictions are not automatically eligible historical screen joins",
        "files": {p.name: sha256(p) for p in args.out.glob("*.csv")}}, indent=2)+"\n")


if __name__ == "__main__":
    main()
