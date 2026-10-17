# Reading guide

Sources checked 17 October 2026. Separate a paper's analyzed release from the
current live database, and distinguish measurements from predictions.

## GenomeSPOT — Barnum et al., 2024

[Paper](https://doi.org/10.1101/2024.03.22.586313) ·
[code and model documentation](https://github.com/cultivarium/GenomeSPOT)

The models infer oxygen tolerance, temperature, salinity and pH from genome/protein
composition. Inputs are organism-level nucleotide and protein FASTAs. Oxygen is
a binary tolerance model; it does not distinguish all respiratory lifestyles.
Continuous outputs include minimum, optimum and maximum growth conditions.

Read the family-level evaluation, compartment-specific features, and random
incompleteness experiment before interpreting predictions. The implementation's
oxygen `error` field is the probability of the reported class. Convert a
not-tolerant class probability q to tolerance probability 1−q. Continuous
`error` is a local cross-validation RMSE estimate, not a calibrated individual CI.
The `is_novel` flag concerns the model's training distribution.

Learning check: reproduce the shipped example, explain the oxygen conversion,
identify missing compartments, and explain why a prediction is not an observed
environmental measurement. The local pilot reproduces the shipped example.

## AllTheBacteria — Hunt et al., version 8, 7 July 2026

[Preprint](https://www.biorxiv.org/content/10.1101/2024.03.08.584059v8) ·
[official documentation](https://allthebacteria.org/docs/)

This community resource standardizes public isolate-genome data and supports
multiple biological applications, including peptide discovery. The analyzed
release contains 2,440,377 high-quality genomes from 11,273 species; the current
resource can be larger. It is not a representative sample of soil organisms.

For this project, follow assembly/sample IDs, release-specific quality metrics,
taxonomy, and annotation provenance. A source genome is the analysis unit;
numerous related isolates and repeated peptide hits are not independent samples.
Start with a stratified subset and publish the selection manifest. Individual
assemblies can be obtained without downloading the entire archive.

## Archaeome — Torres, Wan and de la Fuente-Nunez, 2025

[Nature Microbiology paper](https://www.nature.com/articles/s41564-025-02061-0)

The study mined 18,677 nonredundant reviewed proteins from 233 archaeal organisms
for encrypted antimicrobial peptides, then experimentally tested a selected
set. Reviewed sequences and peptide fragments are not interchangeable with
complete genome-derived proteomes.

For GenomeSPOT, resolve the source organism to a versioned genome, obtain its
corresponding full proteome, and preserve ambiguities. The downloaded
Supplementary Data 2 is a selected hit table, not an exhaustive list of all
233 inputs. Its rows may name more than one organism. The current UniProt query
used in our pilot is not the paper's historical input set.

Predicted source-organism physiology does not establish a peptide's natural
ecological role or its activity under those conditions. Those require separate
evidence. This project does not synthesize, optimize or assay peptides.
