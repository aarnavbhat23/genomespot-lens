"""Generate an explicitly synthetic fixture for notebook execution and tests.

Identifiers, sequences' surrogate IDs, traits, scores, AMR and tree are simulated.
No output from this fixture is evidence about microbial biology.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def make_demo(out, seed=23):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(seed)
    genomes, screens, hits, spot, amr_status, amr_hits = [], [], [], [], [], []
    family_trees = []
    for family in range(12):
        lineage = rng.normal(0, .35)
        species_trees = []
        for species in range(4):
            species_effect = rng.normal(0, .25)
            tips = []
            for isolate in range(8):
                i = len(genomes)
                genome = f"SIM_{i:04d}"
                tips.append(f"{genome}:0.15")
                protein_count = int(rng.integers(1500, 6000))
                proteome_aa = protein_count * rng.uniform(260, 340)
                gc = float(np.clip(.48 + lineage * .08 + rng.normal(0, .04), .25, .75))
                temperature = 35 + 13 * lineage + rng.normal(0, 9)
                ph = float(np.clip(7 + lineage * .3 + rng.normal(0, .6), 3, 10))
                salinity = float(np.clip(1.8 + lineage + rng.normal(0, .8), .1, 10))
                tolerant = rng.random() > .35
                oxygen = float(rng.uniform(.81, .98)) if tolerant else float(rng.uniform(.02, .19))
                n_amr = int(rng.poisson(np.exp(1.2 + species_effect + .3 * lineage)))
                n_screened = int(proteome_aa * rng.uniform(.04, .09))
                eta = 2 + .018 * (temperature - 35) + .15 * (ph - 7) + .09 * n_amr + lineage + species_effect
                mean = n_screened / 100000 * np.exp(eta)
                n_hits = int(rng.negative_binomial(3, 3 / (3 + mean)))
                if i % 29 == 0:
                    n_hits = 0
                genomes.append({"genome_id": genome, "assembly_accession": f"SIM_ASSEMBLY_{i:04d}.1", "cohort": "atb", "domain": "Bacteria", "phylum": f"SIM_Phylum_{family//4}", "family": f"SIM_Family_{family}", "genus": f"SIM_Genus_{family}_{species//2}", "species": f"SIM_Species_{family}_{species}", "protein_count": protein_count, "proteome_aa": proteome_aa, "gc_fraction": gc, "completeness": 91 if i % 37 == 0 else 99, "contamination": 1, "habitat": ["host", "water", "soil"][i % 3]})
                screens.append({"genome_id": genome, "status": "FAIL" if i % 53 == 0 else "COMPLETE", "hits_complete": i % 53 != 0, "n_screened": n_screened, "screen_protocol": "SIMULATED_equal_rule", "exposure_unit": "unique_peptide_sequences", "score_model": "SIMULATED_not_APEX", "score_version": "fixture1", "target_panel_id": "SIM_two_targets"})
                for name, value, error in [("temperature_optimum", temperature, 6), ("ph_optimum", ph, .9), ("salinity_optimum", salinity, 1.7), ("oxygen", "tolerant" if tolerant else "not tolerant", oxygen if tolerant else 1 - oxygen)]:
                    spot.append({"genome_id": genome, "target": name, "value": value, "error": error, "is_novel": i % 41 == 0, "warning": None})
                amr_status.append({"genome_id": genome, "status": "FAIL" if i % 43 == 0 else "PASS"})
                for gene in range(n_amr):
                    amr_hits.append({"genome_id": genome, "element_type": "AMR", "gene_symbol": f"SIM_gene_{gene}", "drug_class": f"SIM_class_{gene % 4}", "locus_id": f"SIM_locus_{gene}"})
                for peptide in range(n_hits):
                    score = 10 ** np.clip(np.log10(26) - .025 * n_amr + rng.normal(0, .2), 0, np.log10(65))
                    cluster = f"SIM_cluster_{rng.integers(0, 25)}"
                    for target, factor in [("SIM_target_A", 1), ("SIM_target_B", 1.15)]:
                        hits.append({"genome_id": genome, "peptide_id": f"SIM_peptide_{i}_{peptide}", "parent_protein_id": f"SIM_parent_{peptide}", "target": target, "predicted_mic_um": score * factor, "mapping_status": "exact", "peptide_cluster_id": cluster})
            species_trees.append("(" + ",".join(tips) + "):0.35")
        family_trees.append("(" + ",".join(species_trees) + "):0.8")
    for name, rows in [("genomes", genomes), ("screening", screens), ("hits", hits), ("genomespot", spot), ("amr_status", amr_status), ("amr_hits", amr_hits)]:
        pd.DataFrame(rows).to_csv(out / f"{name}.csv", index=False)
    (out / "phylogeny.nwk").write_text("(" + ",".join(family_trees) + ");\n")
    spec = {"schema_version": 1, "data_kind": "demo", "seed": seed, "threshold_mode": "any", "threshold_comparison": "lt", "threshold_um": 80, "targets": ["SIM_target_A", "SIM_target_B"], "tree_rooted": True, "tree_rooting_method": "known simulated root", "description": "ALL DATA SIMULATED. Software demonstration only. No biological findings.", "sources": {"all_tables": "scripts/make_study_demo.py"}}
    (out / "bundle.json").write_text(json.dumps(spec, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    make_demo(args.out)
