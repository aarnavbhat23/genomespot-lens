"""Retrieve small published score tables and archive a source-availability audit.

Does not download the 1.25-GB AMR table or model weights. Raw tables stay in runs/.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re

import pandas as pd
import requests

from genomespot_lens.core import sha256

BASE = "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41564-025-02061-0/MediaObjects/41564_2025_2061_"
URLS = {
    "archaeome_all_scores.csv": BASE + "MOESM3_ESM.csv",
    "archaeome_selected_sources.csv": BASE + "MOESM4_ESM.csv",
    "mendeley_files.json": "https://data.mendeley.com/public-api/datasets/d8yzgtdrcp/files?folder_id=root&version=3&$start=0&$limit=1000",
    "carbon_metadata.json": "https://huggingface.co/api/models/HuggingFaceBio/Carbon-A-1.2B",
    "amr_latest_files.json": "https://api.osf.io/v2/nodes/7nwrx/files/osfstorage/69f1e3f3d84f7f2bb2fc8fa3/?page%5Bsize%5D=100",
    "atb_file_index.tsv": "https://osf.io/download/r6gcp/",
    "mendeley_s1.fasta": "https://data.mendeley.com/public-files/datasets/d8yzgtdrcp/files/5c763dac-eccb-4dfd-8750-31bb2f617f6a/file_downloaded",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", required=True)
    parser.add_argument("--audit-out", required=True)
    args = parser.parse_args()
    raw = Path(args.raw_dir)
    raw.mkdir(parents=True, exist_ok=False)
    files = {}
    for name, url in URLS.items():
        response = requests.get(url, timeout=60, headers={"Accept": "application/vnd.mendeley-public-dataset.1+json"} if name == "mendeley_files.json" else {})
        response.raise_for_status()
        (raw / name).write_bytes(response.content)
        files[name] = {"url": url, "sha256": sha256(raw / name), "bytes": len(response.content)}
    full = pd.read_csv(raw / "archaeome_all_scores.csv")
    selected = pd.read_csv(raw / "archaeome_selected_sources.csv")
    carbon = json.loads((raw / "carbon_metadata.json").read_text())
    amr = json.loads((raw / "amr_latest_files.json").read_text())
    mendeley = json.loads((raw / "mendeley_files.json").read_text())
    file_index = pd.read_csv(raw / "atb_file_index.tsv", sep="\t")
    index_matches = file_index.astype(str).apply(lambda col: col.str.contains("peptide|apex|encrypted", case=False, regex=True)).any(axis=1)
    old_s1 = (raw / "mendeley_s1.fasta").read_text(encoding="utf-8-sig")
    report = {
        "checked_utc": datetime.now(timezone.utc).isoformat(), "files": files,
        "archaeome": {"score_rows": len(full), "unique_sequences": int(full.iloc[:, 0].nunique()), "score_columns": list(full.columns), "selected_mapping_rows": len(selected), "selected_mapping_columns": list(selected.columns), "single_taxon_rows": int(selected.ProteinInfo.map(lambda x: len(set(re.findall(r'OX=(\d+)', x))) == 1).sum()), "full_genome_mapping_available_in_score_table": False, "per_genome_screen_denominators_available_in_score_table": False, "scope": "positive score catalog and selected mappings; insufficient for unbiased per-genome richness"},
        "mendeley": {"doi": "10.17632/d8yzgtdrcp.3", "data_s1_content": {"format": "FASTA despite CSV filename", "fasta_records": sum(line.startswith(">") for line in old_s1.splitlines()), "sha256": sha256(raw / "mendeley_s1.fasta")}, "file_inventory": [{"name": x.get("filename", x.get("name")), "size": x.get("content_details", {}).get("size"), "sha256": x.get("content_details", {}).get("sha256_hash"), "url": x.get("content_details", {}).get("download_url")} for x in mendeley]},
        "amr": {"tool_version": "4.2.5", "database_version": "2025-12.03.1", "files": [{"name": x["attributes"]["name"], "bytes": x["attributes"].get("size"), "download": x.get("links", {}).get("download")} for x in amr.get("data", [])], "zero_count_rule": "Only PASS with no AMR determinants is zero; FAIL/NOT DONE/absent remain missing"},
        "carbon_a": {"model": "HuggingFaceBio/Carbon-A-1.2B", "revision": carbon["sha"], "role": "Eukaryotic CDS annotation. Outside the documented scope of bacterial/archaeal gene calling; not an AMP predictor.", "weights_downloaded": False},
        "atb": {"public_genomes_and_amr": True, "master_index": {"url": URLS["atb_file_index.tsv"], "rows": len(file_index), "matched_rows": int(index_matches.sum()), "sha256": sha256(raw / "atb_file_index.tsv"), "searched_terms": ["peptide", "apex", "encrypted"]}, "complete_amp_to_genome_scoring_manifest_located": False, "caveat": "Scoped public search, not a claim that the files do not exist. See study/source-search.md."},
    }
    Path(args.audit_out).write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
