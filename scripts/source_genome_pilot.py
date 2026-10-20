"""Run GenomeSPOT on a bounded set of source organisms from archaeome Data 2.

These are current reference assemblies of reported source taxa, not a
reconstruction of the historical screened proteomes. Counts in Data 2 are
selected-hit counts and must never be used as per-genome AMP richness.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
from importlib.metadata import version
import csv
import gzip
import json
from pathlib import Path
import re
import sys
import time

import requests

from genomespot_lens.cli import clean
from genomespot_lens import UPSTREAM_COMMIT
from genomespot_lens.core import read_fasta, sha256
from genomespot_lens.engine import Engine

SOURCE = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41564-025-02061-0/MediaObjects/41564_2025_2061_MOESM4_ESM.csv"


def download(url, path):
    with requests.get(url, timeout=90, stream=True) as response:
        response.raise_for_status()
        total = 0
        with path.open("wb") as handle:
            for chunk in response.iter_content(1024 * 1024):
                total += len(chunk)
                if total > 50 * 1024 * 1024:
                    raise ValueError("Pilot input exceeds 50 MiB")
                handle.write(chunk)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--limit", type=int, default=6)
    args = parser.parse_args()
    if not 1 <= args.limit <= 20:
        raise ValueError("Pilot limit must be between 1 and 20")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    inputs = out / "inputs"
    inputs.mkdir()
    source = inputs / "archaeome_data2.csv"
    download(SOURCE, source)
    rows = list(csv.DictReader(source.open()))
    counts = Counter()
    for row in rows:
        taxa = set(re.findall(r"OX=(\d+)", row["ProteinInfo"]))
        if len(taxa) == 1:
            counts.update(taxa)
    manifest, prediction_rows = [], []
    run_provenance = {"created_utc": datetime.now(timezone.utc).isoformat(), "genomespot_commit": UPSTREAM_COMMIT,
                      "models": {p.name: sha256(p) for p in Path(args.models).iterdir() if p.is_file()},
                      "software": {p: version(p) for p in ["numpy", "scikit-learn", "hmmlearn", "biopython", "pandas"]}}
    for taxid, n_selected in counts.most_common(args.limit):
        item = {"taxon_id": taxid, "selected_table_rows": n_selected,
                "mapping_scope": "current exact-taxon reference; historical proteome unverified"}
        manifest.append(item)
        try:
            url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxid}/dataset_report?page_size=100"
            response = requests.get(url, timeout=60)
            response.raise_for_status()
            payload = response.json()
            (inputs / f"{taxid}.ncbi.json").write_text(json.dumps(payload, indent=2))
            candidates = {}
            for record in payload.get("reports", []):
                info = record.get("assembly_info", {})
                if str(record.get("organism", {}).get("tax_id")) != taxid or info.get("assembly_level") != "Complete Genome":
                    continue
                accession = record["accession"]
                if not accession.startswith("GCF_"):
                    accession = record.get("paired_accession", "")
                if accession.startswith("GCF_"):
                    candidates[accession] = record
            preferred = {a:r for a,r in candidates.items() if r.get("assembly_info", {}).get("refseq_category") == "reference genome"}
            if len(preferred) == 1:
                candidates = preferred
            if len(candidates) != 1:
                item.update(status="unresolved", reason=f"{len(candidates)} complete exact-taxon RefSeq assemblies")
                continue
            accession, record = next(iter(candidates.items()))
            name = record["assembly_info"]["assembly_name"]
            stem = f"{accession}_{name}"
            digits = accession.split("_")[1].split(".")[0]
            base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{digits[:3]}/{digits[3:6]}/{digits[6:9]}/{stem}"
            fna, faa = inputs / f"{accession}.fna.gz", inputs / f"{accession}.faa.gz"
            download(f"{base}/{stem}_genomic.fna.gz", fna)
            download(f"{base}/{stem}_protein.faa.gz", faa)
            proteins = read_fasta(faa)
            item.update(assembly_accession=accession, organism=record["organism"]["organism_name"], nucleotide_url=f"{base}/{stem}_genomic.fna.gz", protein_url=f"{base}/{stem}_protein.faa.gz", nucleotide_sha256=sha256(fna), protein_sha256=sha256(faa), protein_count=len(proteins), proteome_aa=sum(map(len, proteins.values())))
            engine = Engine(fna, faa, args.models)
            predictions, compartments = engine.predict(sorted(proteins))
            item.update(status="predicted", compartment_counts=compartments)
            for target, prediction in clean(predictions).items():
                prediction_rows.append({"genome_id": accession, "target": target, **prediction})
            print(accession, item["organism"], "predicted", flush=True)
        except (requests.RequestException, ValueError, KeyError) as exc:
            item.update(status="failed", reason=str(exc))
            print(taxid, "failed:", str(exc), file=sys.stderr, flush=True)
        finally:
            (out / "manifest.json").write_text(json.dumps(clean({**run_provenance, "source_url": SOURCE, "source_sha256": sha256(source), "genomes": manifest}), indent=2) + "\n")
            if prediction_rows:
                with (out / "genomespot.csv").open("w") as handle:
                    writer = csv.DictWriter(handle, fieldnames=list(prediction_rows[0]), lineterminator="\n")
                    writer.writeheader()
                    writer.writerows(prediction_rows)
            time.sleep(.4)


if __name__ == "__main__":
    main()
