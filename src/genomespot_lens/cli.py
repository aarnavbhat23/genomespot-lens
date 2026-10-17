"""Command line interface. All output directories must be new."""

import argparse
import csv
from datetime import datetime, timezone
from importlib.metadata import version
import json
from pathlib import Path
import random
import sys

import numpy as np

from . import UPSTREAM_COMMIT, __version__
from .core import map_sequences, matched_sample, numeric_prediction, quality_flags, read_fasta, sha256, validate_ids


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_json(path, value):
    Path(path).write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")


def summarize(records, targets):
    base = records[0]["predictions"]
    selected = records[1]["predictions"]
    summaries = []
    for target in targets:
        baseline = numeric_prediction(target, base[target])
        observed = numeric_prediction(target, selected[target])
        delta = observed - baseline if observed is not None and baseline is not None else None
        for control in ("count_matched", "length_matched"):
            controls = [record["predictions"][target] for record in records if record["kind"] == control]
            shifts = [numeric_prediction(target, p) - baseline for p in controls if baseline is not None and numeric_prediction(target, p) is not None]
            # Descriptive envelope only: not an accuracy CI or hypothesis-test p-value.
            envelope = np.quantile(shifts, [0.025, 0.5, 0.975]).tolist() if shifts else [None] * 3
            summaries.append({
                "target": target, "control": control, "baseline": baseline,
                "selected": observed, "delta": delta,
                "baseline_flags": quality_flags(target, base[target]),
                "selected_flags": quality_flags(target, selected[target]),
                "control_total": len(controls), "control_finite": len(shifts),
                "control_unflagged": sum(not quality_flags(target, p) for p in controls),
                "control_delta_q025": envelope[0], "control_delta_median": envelope[1], "control_delta_q975": envelope[2],
                "fraction_controls_abs_delta_ge_selected": sum(abs(x) >= abs(delta) for x in shifts) / len(shifts) if shifts and delta is not None else None,
            })
    return summaries


def report(out, records, summaries, label):
    fmt = lambda x: "missing" if x is None else f"{x:.4g}"
    lines = [f"# GenomeSPOT Lens — {label}", "", "Exploratory model-input sensitivity audit. Predictions are not measured physiology; subset effects are not causal gene effects.", "", f"Full proteome: {records[0]['protein_count']} proteins. Selected subset: {records[1]['protein_count']} proteins.", "", "| Target | Full input | Selected input | Change | Full / selected flags |", "|---|---:|---:|---:|---|"]
    for s in summaries[::2]:
        flags = " / ".join(", ".join(s[k]) or "none" for k in ("baseline_flags", "selected_flags"))
        lines.append(f"| {s['target']} | {fmt(s['baseline'])} | {fmt(s['selected'])} | {fmt(s['delta'])} | {flags} |")
    lines += ["", "Oxygen numbers are P(oxygen tolerant), not oxygen concentration. Temperature is °C. Other continuous targets retain upstream units.", "", "| Target | Control | Finite / attempted | Unflagged | Control change 2.5–97.5% |", "|---|---|---:|---:|---|"]
    for s in summaries:
        lines.append(f"| {s['target']} | {s['control']} | {s['control_finite']}/{s['control_total']} | {s['control_unflagged']} | {fmt(s['control_delta_q025'])} to {fmt(s['control_delta_q975'])} |")
    lines += ["", "The control envelopes summarize finite predictions, including flagged ones. Missing predictions are counted and excluded from numerical summaries; consult raw records before interpreting results. Envelopes describe random subset variability, not confidence intervals for physiology. Matching uses protein count and fixed length bins; it does not control all biological confounders. Both controls draw from the whole proteome, including the selected group.", "", "DNA is held constant: this is an input-selection audit, not a simulated genome deletion. A full-input prediction is a comparison baseline, not ground truth. Upstream errors for continuous traits are local cross-validation RMSE estimates, not per-genome confidence intervals. Novelty flags indicate unusual model inputs, not scientific discovery.", "", "See provenance.json for checksums and software versions; predictions.json for every subset, protein ID, prediction and warning; summary.tsv for numeric comparisons.", ""]
    (out / "report.md").write_text("\n".join(lines))


def audit(args):
    from .engine import Engine

    proteins = read_fasta(args.proteins)
    read_fasta(args.contigs, protein=False)
    ids = validate_ids(args.subset_ids, proteins)
    if len(ids) == len(proteins):
        raise ValueError("Select fewer proteins than the full proteome")
    if args.replicates < 1 or args.length_bin < 1:
        raise ValueError("Replicates and length-bin must be positive")
    out = Path(args.out)
    if out.exists():
        raise ValueError(f"Output already exists: {out}")
    engine = Engine(args.contigs, args.proteins, args.models)
    available = set(engine.predictor.predict_from_genome(engine.features(sorted(proteins)), engine.models))
    targets = args.targets.split(",")
    if len(set(targets)) != len(targets) or set(targets) - available:
        raise ValueError(f"Targets must be unique and drawn from: {sorted(available)}")
    out.mkdir(parents=True)
    lengths = {k: len(v) for k, v in proteins.items()}
    rng = random.Random(args.seed)
    runs = [("full", sorted(proteins)), ("selected", ids)]
    for name, width in [("count_matched", None), ("length_matched", args.length_bin)]:
        for _ in range(args.replicates):
            runs.append((name, matched_sample(lengths, ids, rng, width)))
    records = []
    for i, (kind, selected) in enumerate(runs):
        predictions, compartments = engine.predict(selected)
        records.append(clean({"run": i, "kind": kind, "protein_count": len(selected), "residue_count": sum(lengths[k] for k in selected), "compartment_counts": compartments, "protein_ids": selected, "predictions": predictions}))
        if i % 20 == 0:
            print(f"{i + 1}/{len(runs)} predictions complete", file=sys.stderr)
    summaries = summarize(records, targets)
    write_json(out / "predictions.json", records)
    write_json(out / "summary.json", summaries)
    with (out / "summary.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: json.dumps(v) if isinstance(v, list) else v for k, v in s.items()} for s in summaries)
    provenance = {
        "lens_version": __version__, "genomespot_commit": UPSTREAM_COMMIT,
        "created_utc": datetime.now(timezone.utc).isoformat(), "python": sys.version,
        "parameters": {k: v for k, v in vars(args).items() if k != "func"},
        "inputs": {k: {"name": Path(v).name, "sha256": sha256(v)} for k, v in {"contigs": args.contigs, "proteins": args.proteins, "subset_ids": args.subset_ids}.items()},
        "models": {p.name: sha256(p) for p in sorted(Path(args.models).iterdir()) if p.is_file()},
        "packages": {k: version(k) for k in ["numpy", "scipy", "scikit-learn", "hmmlearn", "biopython", "pandas", "joblib"]},
        "interpretation": "Descriptive model sensitivity only; no causal inference or measured accuracy claim.",
    }
    write_json(out / "provenance.json", provenance)
    report(out, records, summaries, args.label)
    print(out / "report.md")


def map_subset(args):
    ids, details = map_sequences(read_fasta(args.proteins), read_fasta(args.subset_fasta))
    out = Path(args.out)
    if out.exists():
        raise ValueError(f"Output already exists: {out}")
    out.mkdir(parents=True)
    (out / "subset_ids.txt").write_text("\n".join(ids) + ("\n" if ids else ""))
    write_json(out / "mapping.json", {"full_sha256": sha256(args.proteins), "subset_sha256": sha256(args.subset_fasta), "source_records": len(details), "unique_matched_genes": len(ids), "records": details})
    print(f"{len(ids)} distinct genes matched; inspect {out / 'mapping.json'} for unresolved records")
    if not ids:
        raise ValueError("No unique exact matches; inspect the mapping report")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True)
    mapping = sub.add_parser("map-subset", help="Map a curated protein FASTA to full-proteome IDs by exact sequence")
    mapping.add_argument("--proteins", required=True)
    mapping.add_argument("--subset-fasta", required=True)
    mapping.add_argument("--out", required=True)
    mapping.set_defaults(func=map_subset)
    run = sub.add_parser("audit", help="Compare full/selected proteomes and two matched random controls")
    for name in ("contigs", "proteins", "subset-ids", "models", "out"):
        run.add_argument("--" + name, required=True)
    run.add_argument("--label", default="protein-selection audit")
    run.add_argument("--replicates", type=int, default=100)
    run.add_argument("--seed", type=int, default=23)
    run.add_argument("--length-bin", type=int, default=100)
    run.add_argument("--targets", default="oxygen,temperature_optimum")
    run.set_defaults(func=audit)
    args = parser.parse_args()
    try:
        args.func(args)
    except (ValueError, FileNotFoundError, ImportError) as exc:
        parser.exit(2, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
