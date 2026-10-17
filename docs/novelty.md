# What might be new, and what is already known

This is a scoping review, not proof of priority. No publication claim is made.

| Existing work | Consequence for this project |
|---|---|
| [GenomeSPOT](https://github.com/cultivarium/GenomeSPOT), including `notebooks/analyze_genome_completeness.ipynb` | Random subsampling and genome-wide physiology prediction already exist. Repeating them alone is not the proposed novelty. |
| [metaTraits and porTraits](https://academic.oup.com/nar/article/54/D1/D835/8343513) | Broad trait databases and multi-predictor annotation workflows already exist. A new wrapper or atlas alone is a weak novelty claim. |
| [GenomeSPOT-based growth-rate work](https://www.nature.com/articles/s41467-025-59558-9) | Connecting predicted temperature to other genomic traits already has substantial prior art. |
| [Archaeome study](https://www.nature.com/articles/s41564-025-02061-0) | Selected source proteins create a practical provenance problem worth testing before organism-level inference. |

Candidate contribution: an explicitly tested audit of **nonrandom protein
selection** using protein-count and length-matched controls, mapping attrition,
compartment coverage and model warning rates. A useful result would identify
when annotation or curated-selection practices alter inference beyond random
sampling, and give researchers a reproducible way to diagnose that problem.

Targeted searches included GenomeSPOT with “selection bias”, “protein subsets”,
“annotation bias”, “Swiss-Prot”, “core proteome”, “completeness” and “robustness”.
These searches did not establish that no one has addressed this question.
Before a manuscript, expand to related predictors, citing/cited papers,
proteome representation bias and gene-calling sensitivity. Record search dates,
queries, inclusion criteria and overlapping methods.

Go/no-go: if the cohort shows no systematic effect, publish a useful negative
benchmark and audit tool. Do not manufacture a discovery from two examples,
an outlier flag, or the largest shift found after many unplanned comparisons.
