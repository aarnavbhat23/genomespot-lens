# Current source-genome physiology batch

These are **actual GenomeSPOT predictions**, separate from the simulated
association demonstrations in the three notebooks.

All 233 reconstructed source taxa were checked without using AMP outcomes.
There are 138 unique current complete RefSeq matches, 94 unresolved taxa, and
one bounded query failure for the broad uncultured-archaeon taxon. All 138 input
pairs downloaded successfully and produced 1,380 trait predictions.

- `source_genome_mappings.csv` and `resolution_audit.json`: all taxon statuses,
  versioned current assemblies, source URLs and available NCBI CheckM QC.
- `source_protein_identity.csv`: historical-to-current protein identity/containment,
  full-proteome size and input hashes. All historical source proteins are contained
  in 52 current proteomes; 47 match as full proteins.
- `genomespot.csv`, `genomespot_status.csv`: all predictions/warnings and run statuses.
- `prediction_records/`: per-genome features, predictions and input/model/software
  signatures. The pinned upstream models are referenced by hash, not redistributed.
- `join_checks.csv`, `prediction_audit.json`: identity, QC and primary-trait checks.
  Of 138 genomes, 103 pass QC and 89 have all four primary traits unflagged; the
  intersection with source-protein containment leaves **24 preliminary candidates**.

These 24 candidates are a feasibility result, not a finalized association cohort.
Final analyses require taxonomy, independent habitat, a rooted archaeal tree,
coverage/composition controls and adequate power. No new AMP scoring or measured
activity assay was performed. Current assembly identity is not inferred solely
from a shared organism name.

[Reproduction commands](../../../study/README.md#resolve-and-run-the-reconstructed-source-taxa).
Source genomes/proteomes and full NCBI API responses remain local in `runs/`;
public download URLs and hashes are recorded here.
