"""Download a small, public pilot and run two audits. No bulk collection download.

Run from the repository root after installing .[pilot] and pinned GenomeSPOT.
The selected ATB sample is a download smoke test from the official docs, not
a random or quality-screened research cohort. The archaeal run is a current
reviewed-protein mapping audit, NOT a reconstruction of the paper's inputs.
"""

import argparse
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path
import subprocess
import sys
from urllib.request import urlopen

from genomespot_lens.core import read_fasta, sha256

SOURCES = {
    "atb.fna.gz": "https://allthebacteria-assemblies.s3.eu-west-2.amazonaws.com/SAMD00000344.fa.gz",
    "archaea.fna.gz": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/665/GCF_000091665.1_ASM9166v1/GCF_000091665.1_ASM9166v1_genomic.fna.gz",
    "archaea.faa.gz": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/665/GCF_000091665.1_ASM9166v1/GCF_000091665.1_ASM9166v1_protein.faa.gz",
    "reviewed.faa": "https://rest.uniprot.org/uniprotkb/stream?query=%28organism_id%3A243232%29%20AND%20%28reviewed%3Atrue%29&format=fasta",
    "archaea_assembly.json": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/243232/dataset_report?page_size=5",
    "archaeome_supplementary_data_2.csv": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41564-025-02061-0/MediaObjects/41564_2025_2061_MOESM4_ESM.csv",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models", required=True)
    parser.add_argument("--out", default="runs/pilot")
    parser.add_argument("--replicates", type=int, default=100)
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    data = out / "inputs"
    data.mkdir()
    manifest = []
    for name, url in SOURCES.items():
        path = data / name
        print(f"Downloading {name}", flush=True)
        with urlopen(url, timeout=90) as response, path.open("wb") as handle:
            total = 0
            while chunk := response.read(1024 * 1024):
                total += len(chunk)
                if total > 50 * 1024 * 1024:
                    raise ValueError("Pilot download exceeded 50 MiB per input")
                handle.write(chunk)
        manifest.append({"name": name, "url": url, "bytes": total, "sha256": sha256(path), "retrieved_utc": datetime.now(timezone.utc).isoformat()})
    import pyrodigal
    contigs = read_fasta(data / "atb.fna.gz", protein=False)
    finder = pyrodigal.GeneFinder(meta=False)
    finder.train(*(s.encode() for s in contigs.values()))
    with (data / "atb.faa").open("w") as handle:
        for index, sequence in enumerate(contigs.values()):
            genes = finder.find_genes(sequence.encode())
            for j, gene in enumerate(genes):
                handle.write(f">contig{index}_gene{j}\n{gene.translate().rstrip('*')}\n")
    manifest.append({"name": "atb.faa", "generated_by": f"pyrodigal {pyrodigal.__version__}", "mode": "single genome trained on all contigs, default translation table", "sha256": sha256(data / "atb.faa")})
    (out / "sources.json").write_text(json.dumps(manifest, indent=2) + "\n")
    proteins = read_fasta(data / "atb.faa")
    # Length-selected stress test deliberately creates a biased input. No functional claim.
    selected = sorted(proteins, key=lambda k: (-len(proteins[k]), k))[:max(1, len(proteins) // 5)]
    (data / "atb_longest20_ids.txt").write_text("\n".join(sorted(selected)) + "\n")
    cli = [sys.executable, "-m", "genomespot_lens.cli"]
    subprocess.run(cli + ["map-subset", "--proteins", str(data / "archaea.faa.gz"), "--subset-fasta", str(data / "reviewed.faa"), "--out", str(out / "mapping")], check=True)
    for label, contig, protein, ids in [
        ("ATB_SAMD00000344_longest20_stress_test", data / "atb.fna.gz", data / "atb.faa", data / "atb_longest20_ids.txt"),
        ("M_jannaschii_current_reviewed_exact_matches", data / "archaea.fna.gz", data / "archaea.faa.gz", out / "mapping/subset_ids.txt"),
    ]:
        subprocess.run(cli + ["audit", "--contigs", str(contig), "--proteins", str(protein), "--subset-ids", str(ids), "--models", args.models, "--out", str(out / label), "--label", label, "--replicates", str(args.replicates), "--targets", "oxygen,temperature_optimum,ph_optimum,salinity_optimum"], check=True)


if __name__ == "__main__":
    main()
