"""Run full-proteome GenomeSPOT across a versioned source-genome manifest.

CSV columns: genome_id,contigs,proteins. Paths are relative to the manifest.
Matching checksums resume completed jobs; failures remain explicit missing data.
"""
import argparse
import csv
from datetime import datetime, timezone
from importlib.metadata import version
import json
from pathlib import Path
import re

from genomespot_lens import UPSTREAM_COMMIT
from genomespot_lens.cli import clean
from genomespot_lens.core import read_fasta, sha256
from genomespot_lens.engine import Engine


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--models", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    manifest = Path(args.manifest).resolve()
    records = list(csv.DictReader(manifest.open()))
    if not records or not {"genome_id", "contigs", "proteins"} <= records[0].keys():
        raise ValueError("Manifest needs genome_id,contigs,proteins")
    ids = [row["genome_id"] for row in records]
    if len(ids) != len(set(ids)) or any(not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", x) for x in ids):
        raise ValueError("Genome IDs must be unique, nonempty, safe accession-like strings")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    model_hashes = {p.name: sha256(p) for p in sorted(Path(args.models).iterdir()) if p.is_file()}
    packages = {p: version(p) for p in ["numpy", "scikit-learn", "hmmlearn", "biopython", "pandas"]}
    status, combined = [], []
    for row in records:
        genome_id = row["genome_id"]
        try:
            contigs = (manifest.parent / row["contigs"]).resolve()
            proteins = (manifest.parent / row["proteins"]).resolve()
            signature = {"genomespot_commit": UPSTREAM_COMMIT, "contigs_sha256": sha256(contigs), "proteins_sha256": sha256(proteins), "models": model_hashes, "packages": packages}
            path = out / f"{genome_id}.json"
            if path.exists():
                cached = json.loads(path.read_text())
                if cached.get("signature") != signature:
                    raise ValueError("Cached inputs/software/models differ; choose a new output directory")
                predictions = cached["predictions"]
                state = "RESUMED"
            else:
                protein_records = read_fasta(proteins)
                read_fasta(contigs, protein=False)
                engine = Engine(contigs, proteins, args.models)
                features = engine.features(sorted(protein_records))
                predictions = clean(engine.predictor.predict_from_genome(features, engine.models))
                cached = {"signature": signature, "genome_id": genome_id, "created_utc": datetime.now(timezone.utc).isoformat(), "features": clean(features), "predictions": predictions}
                path.write_text(json.dumps(cached, indent=2, allow_nan=False) + "\n")
                state = "PASS"
            combined.extend({"genome_id": genome_id, "target": target, **p} for target, p in predictions.items())
            status.append({"genome_id": genome_id, "status": state, "reason": ""})
            print(genome_id, state, flush=True)
        except (ValueError, FileNotFoundError, KeyError) as exc:
            status.append({"genome_id": genome_id, "status": "FAIL", "reason": str(exc)})
            print(genome_id, "FAIL", str(exc), flush=True)
    for name, rows, fields in [("genomespot.csv", combined, ["genome_id", "target", "value", "error", "is_novel", "warning", "units"]), ("genomespot_status.csv", status, ["genome_id", "status", "reason"])]:
        # Always replace the aggregate, including with headers only after total
        # failure; an older successful table must never survive as a stale result.
        with (out / name).open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    if any(row["status"] == "FAIL" for row in status):
        raise SystemExit("Some genomes failed; inspect genomespot_status.csv. Failed rows are not zero predictions.")


if __name__ == "__main__":
    main()
