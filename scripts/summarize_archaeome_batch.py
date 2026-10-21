"""Export real source-batch results, preserving identity/QC/trait limitations."""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil

from genomespot_lens.core import quality_flags, sha256

PRIMARY = ["temperature_optimum", "ph_optimum", "salinity_optimum", "oxygen"]


def read(path):
    return list(csv.DictReader(path.open()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    mappings = {r["assembly_accession"]: r for r in read(args.mapping)
                if r["status"] == "resolved_current_exact_taxon"}
    identities = {r["assembly_accession"]: r for r in read(args.inputs / "source_protein_identity.csv")}
    status = read(args.predictions / "genomespot_status.csv")
    if set(mappings) != set(identities) or set(mappings) != {r["genome_id"] for r in status}:
        raise ValueError("Incomplete batch: retain all resolved mapping/input/prediction statuses")
    args.out.mkdir(parents=True, exist_ok=True)
    raw = args.out / "prediction_records";raw.mkdir(exist_ok=True)
    join_checks, flag_counts, prediction_records = [], Counter(), {}
    for row in status:
        accession = row["genome_id"]
        source, identity = mappings[accession], identities[accession]
        try:
            qc_pass = float(source["completeness"]) >= 95 and float(source["contamination"]) <= 5
        except ValueError:
            qc_pass = False
        predicted = row["status"] in {"PASS", "RESUMED"}
        flags = {}
        if predicted:
            path = args.predictions / f"{accession}.json"
            record = json.loads(path.read_text())
            for target in PRIMARY:
                flags[target] = quality_flags(target, record["predictions"].get(target, {}))
                flag_counts.update(f"{target}:{flag}" for flag in flags[target])
            shutil.copy2(path, raw / path.name)
            prediction_records[path.name] = sha256(path)
        all_source = identity["all_historical_sources_contained"] == "True"
        traits_pass = predicted and not any(flags.values())
        join_checks.append({"taxon_id": source["taxon_id"], "assembly_accession": accession,
            "historical_sources_contained": all_source,
            "historical_sources_exact_full_proteins": identity["exact_whole_protein_matches"] == identity["historical_distinct_source_proteins"],
            "genome_qc_pass": qc_pass, "prediction_success": predicted,
            "all_four_primary_traits_unflagged": traits_pass,
            "preliminary_historical_join_candidate": all_source and qc_pass and traits_pass,
            "trait_flags": json.dumps(flags, sort_keys=True)})
    with (args.out / "join_checks.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(join_checks[0]), lineterminator="\n")
        writer.writeheader();writer.writerows(join_checks)
    for name in ["genomespot.csv", "genomespot_status.csv"]:
        shutil.copy2(args.predictions / name, args.out / name)
    for name in ["source_protein_identity.csv", "download_status.csv"]:
        shutil.copy2(args.inputs / name, args.out / name)
    shutil.copy2(args.inputs / "audit.json", args.out / "download_audit.json")
    counts = {key: sum(r[key] for r in join_checks) for key in [
        "historical_sources_contained", "historical_sources_exact_full_proteins", "genome_qc_pass",
        "prediction_success", "all_four_primary_traits_unflagged", "preliminary_historical_join_candidate"]}
    audit = {"created_utc": datetime.now(timezone.utc).isoformat(), "source_taxa_considered": 233,
        "current_resolved_genomes": len(mappings), "prediction_rows": len(read(args.out / "genomespot.csv")),
        "counts": counts, "trait_flag_counts": dict(flag_counts),
        "scope": "Actual whole-genome/full-proteome GenomeSPOT predictions. Current-reference technical batch, not a completed historical association analysis.",
        "candidate_caveat": "Preliminary join candidates still require documented taxonomy, independent habitat, a rooted archaeal tree, composition/coverage controls and a frozen analysis. Failed QC and flagged traits are retained, not set to zero.",
        "selection": "All 233 reconstructed taxa checked; no AMP outcome used to choose source assemblies",
        "mapping_sha256": sha256(args.mapping),
        "files": {p.name: sha256(p) for p in sorted(args.out.glob("*.csv"))},
        "prediction_records": prediction_records}
    (args.out / "prediction_audit.json").write_text(json.dumps(audit, indent=2)+"\n")
    print(json.dumps({"resolved_genomes": len(mappings), "prediction_rows": audit["prediction_rows"], **counts}, indent=2))


if __name__ == "__main__":
    main()
