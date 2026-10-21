"""Prepare a REAL temperature/AMP study from verified historical protein subsets.

No regressions or candidate selection occur here. TEMPURA and GTDB are separate
source resources with their own licenses; raw inputs remain under runs/.
"""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import html
import json
from pathlib import Path
import re
import unicodedata
from datetime import datetime, timezone

from Bio import Phylo

from genomespot_lens.archaeome_archive import distinct_fragment_count
from genomespot_lens.core import quality_flags, read_fasta, sha256

AA = "ACDEFGHIKLMNPQRSTVWY"


def rows(path):
    return list(csv.DictReader(path.open()))


def write(path, records):
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator="\n")
        w.writeheader();w.writerows(records)


def parse_tempura(path):
    raw = path.read_text().strip()
    wrapped = re.fullmatch(r"/\*\*/callback\((.*)\);?", raw, re.S)
    data = json.loads(wrapped.group(1) if wrapped else raw)
    cols = ["entry", "species", "taxid", "strain", "domain", "phylum", "class", "order", "family", "genus", "assembly", "gc", "size", "rrna", "rrna_gc", "tmin", "topt", "topt_low", "topt_high", "tmax", "trange"]
    result = []
    for record in data["rows"]:
        if len(record["cell"]) != len(cols):
            raise ValueError("TEMPURA schema changed")
        row = {k: html.unescape(re.sub("<[^>]+>", "", str(v))) for k,v in zip(cols, record["cell"])}
        row["entry_id"] = record["id"]
        result.append(row)
    if len(result) != int(data["total"]):
        raise ValueError("Incomplete TEMPURA page")
    return result


def normalized_strain(value):
    value = str(value).replace("∆", "Delta ").replace("Δ", "Delta ")
    return re.sub("[^a-z0-9]", "", unicodedata.normalize("NFKD", value).lower())


def strain_supported(temperature_rows, source_organism, ncbi, verified_aliases=()):
    """Corroborate culture labels; an assembly field alone can contain database errors."""
    original = re.search(r"\(strain ([^)]+)\)", source_organism)
    aliases = original.group(1).split("/") if original else []
    aliases += ncbi.get("assembly_info", {}).get("biosample", {}).get("strain", "").split(";")
    aliases += [ncbi.get("organism", {}).get("infraspecific_names", {}).get("strain", "")]
    aliases += list(verified_aliases)
    aliases = {normalized_strain(x) for x in aliases if normalized_strain(x)}
    return bool(temperature_rows) and all(normalized_strain(t["strain"]) in aliases for t in temperature_rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--raw",type=Path,required=True)
    p.add_argument("--reconstruction",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    root=Path(__file__).resolve().parents[1]
    batch=root/"examples/amp-study/reconstructed-genomes"
    verified_aliases=json.loads((root/"study/temperature/strain-aliases.json").read_text())
    source=rows(a.reconstruction/"protein_registry.csv")
    fasta=read_fasta(a.reconstruction/"reviewed_proteins.fasta")
    by_tax=defaultdict(dict)
    for r in source:
        seq=fasta[r["accession"]]
        if hashlib.sha256(seq.encode()).hexdigest()!=r["sequence_sha256"]:
            raise ValueError("Historical FASTA/registry mismatch")
        by_tax[r["taxon_id"]].setdefault(r["sequence_sha256"],r)
    mapped=defaultdict(set)
    for r in rows(a.reconstruction/"positive_parent_mappings.csv"):
        mapped[(r["taxid"],r["accession"])].add(r["catalog_row"])
    clusters={}
    for cutoff in [30,50]:
        with (a.raw/f"homology{cutoff}_cluster.tsv").open() as fh:
            clusters[cutoff]={member:rep for rep,member in csv.reader(fh,delimiter="\t")}
        if set(clusters[cutoff]) != set(fasta):
            raise ValueError("Homology cluster membership incomplete")
    temperatures=parse_tempura(a.raw/"tempura.json")
    with gzip.open(a.raw/"gtdb-ar53.tree.gz","rt") as fh:tree=Phylo.read(fh,"newick")
    tips={x.name.removeprefix("RS_").removeprefix("GB_"):x.name for x in tree.get_terminals()}
    genomes,proteins,phenotypes,exclusions=[],[],[],[]
    for r in rows(batch/"source_genome_mappings.csv"):
        if r["status"]!="resolved_current_exact_taxon":
            exclusions.append({"taxon_id":r["taxon_id"],"stage":"assembly","reason":r["reason"]});continue
        tax,asm=r["taxon_id"],r["assembly_accession"]
        ncbi=json.loads((root/f"runs/archaeome-genome-resolution/ncbi/{tax}.selected.json").read_text())
        aliases={asm,ncbi["accession"],ncbi.get("paired_accession","")}-{ "" }
        temp=[t for t in temperatures if t["domain"]=="Archaea" and t["assembly"] in aliases]
        observed={float(t["topt"]) for t in temp if t["topt"]!="-"}
        strain_ok=strain_supported(temp,r["source_organism"],ncbi,verified_aliases.get(asm,{}).get("accepted_TEMPURA_strain_labels",[]))
        tmean=next(iter(observed)) if len(observed)==1 and strain_ok else None
        if temp and not strain_ok:
            exclusions.append({"taxon_id":tax,"stage":"strain","reason":"TEMPURA culture label not corroborated by source/NCBI strain aliases"})
        if len(observed)!=1:
            exclusions.append({"taxon_id":tax,"stage":"temperature","reason":"no exact assembly-pair match" if not observed else "conflicting curated temperatures"})
        current=read_fasta(root/f"runs/archaeome-current-inputs/inputs/{asm}.faa.gz")
        present=set(current.values());valid=[];positive_union=set()
        for record in by_tax[tax].values():
            seq=fasta[record["accession"]]
            if seq not in present:continue
            candidates=mapped[(tax,record["accession"])];positive_union.update(candidates)
            exposure=distinct_fragment_count([seq]);valid.append(seq)
            proteins.append({"taxon_id":tax,"genome_id":asm,"parent_id":record["accession"],
                "sequence_sha256":record["sequence_sha256"],"length_aa":len(seq),
                "n_screened":exposure,"amp_count":len(candidates),
                "homology30":clusters[30][record["accession"]],"homology50":clusters[50][record["accession"]],
                **{f"aa_{aa}":seq.count(aa)/len(seq) for aa in AA}})
        pred=json.loads((batch/"prediction_records"/f"{asm}.json").read_text())["predictions"]["temperature_optimum"]
        try:qc=float(r["completeness"])>=95 and float(r["contamination"])<=5
        except ValueError:qc=False
        concat="".join(valid);full="".join(current.values())
        if not valid:
            exclusions.append({"taxon_id":tax,"stage":"identity","reason":"no exact full historical protein in current annotation"})
        reasons=[]
        if not qc:reasons.append("missing_or_failed_genome_qc")
        if len(valid)<20:reasons.append("fewer_than_20_verified_parents")
        if tmean is None:reasons.append("unresolved_curated_strain_temperature")
        if asm not in tips:reasons.append("no_exact_marker_tree_tip")
        genomes.append({"genome_id":asm,"taxon_id":tax,"organism":r["source_organism"],
            "genome_qc_pass":qc,"matched_source_proteins":len(valid),
            "historical_unique_proteins":len(by_tax[tax]),"proteome_aa":len(full),"proteome_proteins":len(current),
            "reviewed_aa_fraction":len(concat)/len(full),"n_screened":distinct_fragment_count(valid),
            "amp_count":len(positive_union),"observed_temperature":tmean,"temperature_strain_supported":strain_ok,
            "predicted_temperature":pred["value"],"predicted_temperature_unflagged":not quality_flags("temperature_optimum",pred),
            "gtdb_tip":tips.get(asm,""),"tempura_match":"exact_versioned_assembly_pair" if tmean is not None else "unresolved",
            "primary_eligible":qc and len(valid)>=20 and tmean is not None and asm in tips,"primary_exclusion_reasons":";".join(reasons),
            **{f"aa_{aa}":concat.count(aa)/len(concat) if concat else 0 for aa in AA}})
        if temp:
            phenotypes.append({"genome_id":asm,"matched_assembly_ids":";".join(sorted({t['assembly'] for t in temp})),
                "tempura_entry_ids":";".join(str(t["entry_id"]) for t in temp),
                "species":";".join(sorted({t["species"] for t in temp})),"strain":";".join(sorted({t["strain"] for t in temp})),
                "topt":tmean,"strain_corroborated":strain_ok,"reported_topt":";".join(sorted({t["topt"] for t in temp})),"tmin":";".join(sorted({t["tmin"] for t in temp})),"tmax":";".join(sorted({t["tmax"] for t in temp})),
                "topt_low":";".join(sorted({t["topt_low"] for t in temp})),"topt_high":";".join(sorted({t["topt_high"] for t in temp})),
                "source_url":"https://togodb.org/db/tempura","license":"CC BY-NC 3.0"})
        print(asm,'verified parents',len(valid),'TEMPURA',tmean,'primary',genomes[-1]['primary_eligible'],flush=True)
    write(a.out/"genomes.csv",genomes);write(a.out/"proteins.csv",proteins)
    write(a.out/"phenotypes.csv",phenotypes);write(a.out/"exclusions.csv",exclusions)
    # Keep every available exact assembly tip, so prespecified cohort sensitivities
    # share one marker tree. Pruning preserves intervening branch lengths.
    selected={g['gtdb_tip']:g['genome_id'] for g in genomes if g['gtdb_tip']}
    def retain(clade):
        if clade.is_terminal():
            if clade.name not in selected:return None
            clade.name=selected[clade.name];return clade
        clade.clades=[child for c in clade.clades if (child:=retain(c)) is not None]
        if not clade.clades:return None
        if len(clade.clades)==1:
            child=clade.clades[0];child.branch_length=(child.branch_length or 0)+(clade.branch_length or 0);return child
        return clade
    tree.root=retain(tree.root);tree.root.branch_length=0;tree.rooted=True
    Phylo.write(tree,a.out/"phylogeny.nwk","newick",format_branch_length="%1.10f")
    audit={"created_utc":datetime.now(timezone.utc).isoformat(),"data_kind":"real","protocol_commit":"bdb6752",
        "genomes":len(genomes),"verified_parent_rows":len(proteins),"primary_genomes":sum(g['primary_eligible'] for g in genomes),
        "phenotype_assembly_matches":len(phenotypes),"phenotype_strain_supported":sum(g["observed_temperature"] is not None for g in genomes),"tempura_rows":len(temperatures),"tree_exact_genome_matches":len(selected),
        "homology_command":"mmseqs easy-cluster reviewed_proteins.fasta PREFIX TMP --min-seq-id {0.3,0.5} -c 0.8 --cov-mode 0 --cluster-mode 2 --threads 4 -s 7.5",
        "homology_version":"18-8cc5c","tree_source":"https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_r226.tree.gz",
        "tree_rooting":"Published GTDB r226 ar53 root; unary paths collapsed preserving branch lengths; no independent rerooting",
        "tempura_source":"https://togodb.org/togodb/fetch/1376?rp=10000&page=1",
        "licenses":{"TEMPURA":"CC BY-NC 3.0","GTDB_tree":"CC BY-SA 4.0","UniProt":"CC BY 4.0"},
        "verified_strain_aliases_sha256":sha256(root/"study/temperature/strain-aliases.json"),
        "raw_input_sha256":{n:sha256(a.raw/n) for n in ["tempura.json","gtdb-ar53.tree.gz","homology30_cluster.tsv","homology50_cluster.tsv"]},
        "outputs":{x.name:sha256(x) for x in a.out.iterdir()}}
    (a.out/"provenance.json").write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({k:audit[k] for k in ['genomes','verified_parent_rows','primary_genomes','phenotype_assembly_matches','tree_exact_genome_matches']},indent=2))


if __name__=="__main__":main()
