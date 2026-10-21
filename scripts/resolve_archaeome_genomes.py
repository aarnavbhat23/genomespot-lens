"""Resolve ALL recovered archaeome taxa to current exact-taxon complete RefSeq assemblies.

Selection uses taxonomy and assembly metadata, never AMP outcomes. Unresolved,
ambiguous, failed and duplicate mappings remain explicit. This does not establish
historical protein sequence identity; use the retained metadata for verification.
"""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import time

import requests

from genomespot_lens.core import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--taxa", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    raw = args.out / "ncbi"
    raw.mkdir()
    rows, seen = [], set()
    taxa = list(csv.DictReader(args.taxa.open()))
    session = requests.Session()
    for n, taxon in enumerate(taxa, 1):
        taxid = taxon["taxon_id"]
        item = {"taxon_id": taxid, "source_organism": taxon["organism"], "status": "failed",
                "assembly_accession": "", "organism": "", "completeness": "", "contamination": "",
                "contigs_url": "", "proteins_url": "", "reason": ""}
        try:
            candidates, reports, token = {}, [], ""
            for page in range(10):
                url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxid}/dataset_report"
                params = {"page_size": 100}
                if token:
                    params["page_token"] = token
                response = session.get(url, params=params, timeout=60)
                response.raise_for_status()
                payload = response.json()
                reports.extend(payload.get("reports", []))
                token = payload.get("next_page_token", "")
                if not token:
                    break
                time.sleep(.36)
            (raw / f"{taxid}.json").write_text(json.dumps({"reports": reports, "remaining_page_token": token}, indent=2))
            if token:
                raise ValueError("More than 1,000 descendant assemblies; incomplete candidate enumeration")
            for report in reports:
                info = report.get("assembly_info", {})
                if str(report.get("organism", {}).get("tax_id")) != taxid or info.get("assembly_level") != "Complete Genome" or info.get("assembly_status") != "current":
                    continue
                accession = report["accession"]
                if not accession.startswith("GCF_"):
                    accession = report.get("paired_accession", "")
                if accession.startswith("GCF_"):
                    # Prefer the actual RefSeq report if both GenBank/RefSeq are returned.
                    if accession not in candidates or report["accession"].startswith("GCF_"):
                        candidates[accession] = report
            preferred = {a:r for a,r in candidates.items()
                         if r.get("assembly_info", {}).get("refseq_category") == "reference genome"}
            if len(preferred) == 1:
                candidates = preferred
            if len(candidates) != 1:
                item.update(status="unresolved", reason=f"{len(candidates)} eligible exact-taxon complete RefSeq assemblies")
            else:
                accession, report = next(iter(candidates.items()))
                if accession in seen:
                    item.update(status="duplicate", reason=f"Assembly {accession} already assigned to another source taxon")
                else:
                    seen.add(accession)
                    stem = accession + "_" + report["assembly_info"]["assembly_name"]
                    digits = accession.split("_")[1].split(".")[0]
                    base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/{digits[:3]}/{digits[3:6]}/{digits[6:9]}/{stem}/{stem}"
                    qc = report.get("checkm_info", {})
                    item.update(status="resolved_current_exact_taxon", assembly_accession=accession,
                        organism=report["organism"]["organism_name"], completeness=qc.get("completeness", ""),
                        contamination=qc.get("contamination", ""), contigs_url=base+"_genomic.fna.gz",
                        proteins_url=base+"_protein.faa.gz", reason="Current reference; historical source-protein identity still requires checking")
                    (raw / f"{taxid}.selected.json").write_text(json.dumps(report, indent=2))
        except (requests.RequestException, ValueError, KeyError) as exc:
            item.update(status="failed", reason=str(exc))
        rows.append(item)
        with (args.out / "source_genome_mappings.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(item), lineterminator="\n")
            writer.writeheader();writer.writerows(rows)
        print(n, "/", len(taxa), taxid, item["status"], item["assembly_accession"], flush=True)
        time.sleep(.36)
    (args.out / "audit.json").write_text(json.dumps({
        "created_utc": datetime.now(timezone.utc).isoformat(), "taxa_sha256": sha256(args.taxa),
        "input_taxa": len(taxa), "statuses": {s: sum(r["status"]==s for r in rows) for s in sorted({r["status"] for r in rows})},
        "selection": "All recovered taxa; exact current NCBI taxon, current complete RefSeq assembly; unique reference preferred; duplicates excluded",
        "scope": "Current-reference mapping, not verified historical protein identity. No AMP outcome used in selection.",
        "mapping_sha256": sha256(args.out / "source_genome_mappings.csv")}, indent=2)+"\n")


if __name__ == "__main__":
    main()
