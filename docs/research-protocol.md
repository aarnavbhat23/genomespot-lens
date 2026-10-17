# Research protocol and route to the mentor's analysis

This document specifies proposed work. The two-genome pilot is complete;
the cohort studies below have not been run. Suggested sizes are planning
choices, not a power calculation or evidence of sufficient sample size.

## First standalone project: genetics and microbial physiology

Working question: **When does nonrandom representation of a genome's proteins
alter inferred oxygen tolerance and temperature preference beyond matched
random sampling?**

The public benefit is a reusable check for researchers using curated proteins,
functional gene panels, or pangenome partitions with genome-level models.
The first release accepts arbitrary supplied protein groups. Identifying core,
accessory or functional groups requires independent annotation and is not
implemented by Lens itself.

1. Freeze the protocol and broaden the novelty review before looking for
   significant results. Keep oxygen and optimum temperature as primary targets;
   treat pH and salinity as secondary until the ecology study.
2. Assemble an initial cohort of roughly 100 versioned reference genomes, aiming
   for broad bacterial/archaeal and family coverage. Prefer complete genomes;
   if drafts are included, predefine completeness ≥95% and contamination ≤5%
   and record the estimation software. These thresholds are proposed criteria.
3. Recover reviewed-protein membership from a fixed UniProt release and map to
   the same strain and assembly. Preserve exact/unmatched/ambiguous mapping
   rates. Set a high mapping-coverage requirement in advance (e.g. ≥95% of
   source proteins), then perform sensitivity analyses. The archaeal pilot's
   ~66% exact mapping does not satisfy that proposed primary-analysis gate.
4. If pursuing core/accessory genetics, construct independently defined
   pangenomes within species with adequate isolates; never label genes “core”
   based on the GenomeSPOT result. Analyze this as a separate, prespecified
   comparison. Do not combine arbitrary panels into one biological class.
5. Compare full-input predictions, selected groups and both matched controls.
   Use at least 100 controls for initial screening, and increase only where
   precision justifies it. Record selection fraction, residues, length bins,
   protein localization, missing features, warnings and classification flips.
6. Analyze at the genome level. Controls from one genome are Monte Carlo
   replicates, not independent organisms. Use paired effect estimates and
   family-blocked bootstrap intervals; examine taxonomy, quality and subset
   size as confounders. Use a phylogenetic comparative model if a trustworthy
   tree is available. Publish unfiltered and quality-filtered summaries.
7. Evaluate accuracy only against independent strain-resolved experimental
   physiology. Establish training-data overlap and exclude overlap from any
   claimed external-validation result. Without independent measurements,
   conclusions must remain about stability, not correctness.
8. Report primary effects even when small or null. If formal multiple hypothesis
   tests are added, preregister hypotheses and control false discovery rate.
   Replicate across another genome/annotation release and held-out families.

Expected figures: cohort flow and exclusions; selected-versus-full changes
with matched control distributions; warning/missingness rates by compartment
coverage; family-stratified effect estimates; independent validation if available.
Manuscript decision follows cohort evidence and prior-art review, not repository
creation. A methods/resource paper may be realistic; acceptance is not guaranteed.

## Then: AllTheBacteria and archaeome physiology

1. **Resolve cohorts.** Freeze an ATB release and a small taxonomically balanced
   discovery panel, initially 100–500 genomes. For the archaeome study recover
   the original list of all 233 source organisms and protein accessions from
   archived input data or authors. Supplementary Data 2 alone is insufficient.
2. **Resolve genomes.** Store source paper ID, source organism string, strain,
   taxon ID, assembly accession/version, protein annotation release, source URL,
   checksum and mapping confidence. Keep unresolved and one-to-many mappings
   visible. Do not force species-name-only matches.
3. **Quality control.** Obtain assembly/contamination/completeness metadata,
   taxonomy and habitat records. De-duplicate by accession and a stated
   species/strain rule. Keep clinically overrepresented lineages from dominating.
4. **Predict from full inputs.** Run GenomeSPOT on each complete organism input,
   preserve all ten predictions and flags, and use Lens to investigate input
   sensitivity where source-protein selection is relevant. Cache features and
   log failures. Scale only after measuring runtime, disk use and error rates.
5. **Separate physiology and environment.** A model predicts an organism's
   growth preference. Claims about soil or in-situ pH require independent,
   traceable environmental measurements. Metadata absence is missingness, not
   evidence that a sample was non-soil. Habitat labels and measured pH must
   remain distinct fields.
6. **Discovery with controls.** Test prespecified taxon/habitat contrasts with
   phylogenetic and sampling controls. If relating physiology to published
   peptide hit counts, normalize by the number/length of screened parent
   proteins and the search procedure; preserve multi-source peptide mappings.
   Shared amino-acid inputs can create circular associations between models.
   Predictions do not demonstrate expression, secretion or natural activity.
7. **Deliver.** Release the manifest, reproducible workflow, uncertainty-aware
   tables, exclusions, figures and negative findings. Use independent data or
   collaborator-led measurements to test any proposed biological mechanism.

## Practical sequence

Suggested week 1: learn the models, inspect this pilot, finalize the question
and identify exact source releases. Weeks 2–3: implement the cohort manifest
and run the selection-bias study. Weeks 4–6: replicate, validate and decide
whether a manuscript is supported. The ATB/archaeome expansion follows the
standalone study, while accession recovery can proceed alongside it.

Compute and time depend on cohort and hardware. No full ATB download, paid
compute, wet-lab experiment, author outreach or mentor message is scheduled by
this repository. Benchmark a small batch before planning the full collection.
