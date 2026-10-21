"""Build exactly three annotated notebooks. Execution is handled separately."""
from pathlib import Path
import textwrap

import nbformat as nbf

ROOT = Path(__file__).resolve().parents[1]


def md(text):
    return nbf.v4.new_markdown_cell(textwrap.dedent(text).strip())


def code(text):
    return nbf.v4.new_code_cell(textwrap.dedent(text).strip())


SETUP = '''
from pathlib import Path
import os, sys, json, platform, subprocess
from importlib.metadata import version
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, Markdown
get_ipython().run_line_magic("matplotlib", "inline")

ROOT = next(p for p in [Path.cwd(), *Path.cwd().parents]
            if (p / "src/genomespot_lens/amp_study.py").exists())
sys.path.insert(0, str(ROOT / "src"))
from genomespot_lens.amp_study import (
    TRAITS, load_bundle, prepare_cohort, analysis_subset, standardized,
    fit_nb, effect_table, grouped_ecology_cv, tree_covariance,
    build_phylogenetic_nb,
)

# Explicit demonstration mode. Change to "real" only after preparing the contract.
# No automatic fallback to simulated data is permitted.
MODE = os.environ.get("AMP_STUDY_MODE", "demo")
assert MODE in {"demo", "real"}
COHORT = os.environ.get("AMP_STUDY_COHORT", "atb")
DATA_DIR = (ROOT / "examples/amp-study/demo") if MODE == "demo" else Path(
    os.environ.get("AMP_STUDY_DATA", str(ROOT / "data/amp-study")))
RUN_PHYLO_MCMC = os.environ.get("AMP_STUDY_PHYLO_MCMC", "0") == "1"
SEED = 23
LABEL = "SIMULATED DATA — software demonstration" if MODE == "demo" else "REAL COHORT — observational associations"
plt.rcParams.update({"figure.dpi": 120, "axes.spines.top": False,
                     "axes.spines.right": False, "font.size": 10})
display(Markdown(f"**{LABEL}**"))
print("Python", platform.python_version(), "| cohort:", COHORT, "| bundle:", "examples/amp-study/demo" if MODE == "demo" else DATA_DIR.name)
FIG_DIR = ROOT / "runs/amp-study/figures" / MODE
FIG_DIR.mkdir(parents=True, exist_ok=True)
def save_figure(figure, name):
    figure.savefig(FIG_DIR / f"{name}.png", dpi=180, bbox_inches="tight")
    figure.savefig(FIG_DIR / f"{name}.svg", bbox_inches="tight")
'''

AUDIT = '''
# Actual published-source audit, independent of the simulated analysis below.
source_audit = json.loads((ROOT / "examples/amp-study/source-audit.json").read_text())
display(pd.DataFrame([
    {"source": "Archaeome complete positive-score catalog", "records": source_audit["archaeome"]["score_rows"], "genome richness ready": False},
    {"source": "Archaeome selected source mappings", "records": source_audit["archaeome"]["selected_mapping_rows"], "genome richness ready": False},
]))
source_dir = ROOT / "examples/amp-study/source-genomes"
if source_dir.exists():
    # Genuine GenomeSPOT outputs. They are not joined to invented AMP counts.
    real_predictions = pd.read_csv(source_dir / "genomespot.csv")
    display(real_predictions[real_predictions.target.isin(TRAITS.values())])
reconstruction_dir = ROOT / "examples/amp-study/archaeome-reconstruction"
reconstruction = json.loads((reconstruction_dir / "audit.json").read_text())
display(pd.DataFrame({"reported": reconstruction["expected_counts"],
                      "reconstructed": reconstruction["observed_counts"]}))
source_taxa = pd.read_csv(reconstruction_dir / "source_taxa.csv")
display(source_taxa[["reviewed_protein_records", "distinct_screened_fragments",
                     "distinct_high_scoring_peptides"]].describe())
print("Actual reconstructed taxa:", len(source_taxa), "| zero catalog hits:",
      reconstruction["taxa_with_zero_catalog_hits"], "| all published positives map:",
      reconstruction["all_positive_peptides_mapped"])
print("These are source-taxon reviewed-protein screens; use verified assembly joins for genome models.")
batch_dir = ROOT / "examples/amp-study/reconstructed-genomes"
if (batch_dir / "prediction_audit.json").exists():
    batch = json.loads((batch_dir / "prediction_audit.json").read_text())
    display(pd.Series(batch["counts"], name="ACTUAL current-reference batch checks"))
    print("Actual GenomeSPOT prediction rows:", batch["prediction_rows"])
    print(batch["candidate_caveat"])
print("ATB Results and Methods disagree on the MIC selection rule; resolve before a real ATB run.")
print("Source-screen linkage: see study/source-search.md and study/data-contract.md")
'''

LOAD = '''
tables, spec, provenance, tree_path = load_bundle(DATA_DIR, COHORT, expected_kind=MODE)
cohort, candidate_table, hit_audit = prepare_cohort(tables, spec)
display(pd.Series(hit_audit, name="hit-table audit"))
display(pd.DataFrame({
    "genomes": [len(cohort)],
    "completed_screens": [int(cohort.screen_eligible.sum())],
    "quality_pass": [int(cohort.genome_qc_pass.sum())],
    "all_four_traits_unflagged": [int(cohort.all_traits_unflagged.sum())],
    "amr_successful": [int(cohort.amr_status.eq("PASS").sum())],
}))
assert not cohort.loc[~cohort.screen_eligible, "amp_count"].notna().any()
assert not cohort.loc[~cohort.amr_status.eq("PASS"), "amr_gene_count"].notna().any()
'''


def export_cell(name, outputs):
    assignments = "\n".join(f'{value}.to_csv(OUT / "{filename}", index=False)' for filename, value in outputs.items())
    return code(f'''
OUT = ROOT / "runs/amp-study/{name}" / MODE
OUT.mkdir(parents=True, exist_ok=True)
{assignments}
provenance["analysis"] = "{name}"
provenance["seed"] = SEED
provenance["run_phylogenetic_mcmc"] = RUN_PHYLO_MCMC
provenance["git_commit"] = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
provenance["git_dirty"] = bool(subprocess.check_output(["git", "-C", str(ROOT), "status", "--porcelain"], text=True).strip())
provenance["bundle_settings"] = spec
provenance["packages"] = {{p: version(p) for p in ["numpy", "pandas", "statsmodels", "pymc", "arviz", "biopython", "nbformat", "nbclient"]}}
from genomespot_lens.core import sha256
provenance["analysis_code_sha256"] = {{str(p.relative_to(ROOT)): sha256(p) for p in [ROOT / "src/genomespot_lens/amp_study.py", ROOT / "scripts/build_study_notebooks.py"]}}
(OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\\n")
print("Saved results:", OUT.relative_to(ROOT))
display(Markdown("**Demonstration complete; simulated coefficients are not biological findings.**"
                 if MODE == "demo" else "**Review diagnostics, exclusions and sensitivity analyses before interpreting associations.**"))
''')


def phylo_cell(design_formula, analysis_name):
    return code(f'''
import patsy
import pymc as pm
if not spec.get("tree_rooted"):
    raise ValueError("A documented rooted tree with branch lengths is required")
phylo_cov = tree_covariance(tree_path, data.genome_id.tolist())
X_phylo = patsy.dmatrix({design_formula!r}, data, return_type="dataframe", NA_action="raise")
phylo_model = build_phylogenetic_nb(data, X_phylo, phylo_cov)
print("Built phylogenetic negative-binomial model:", len(data), "genomes,", X_phylo.shape[1], "fixed effects")
# A prior-predictive check exercises the actual mixed-model likelihood.
with phylo_model:
    prior_check = pm.sample_prior_predictive(samples=20, random_seed=SEED, var_names=["amp_count"])
print("Prior predictive dimensions:", dict(prior_check.prior_predictive.sizes))
if RUN_PHYLO_MCMC:
    import arviz as az
    with phylo_model:
        phylo_fit = pm.sample(draws=1000, tune=1000, chains=4, cores=2,
                              target_accept=.95, random_seed=SEED,
                              idata_kwargs={{"log_likelihood": True}})
        posterior_check = pm.sample_posterior_predictive(phylo_fit, random_seed=SEED)
    diagnostics = az.summary(phylo_fit, var_names=["beta", "phylo_sd", "nb_shape"])
    display(diagnostics)
    print("Divergences:", int(phylo_fit.sample_stats.diverging.sum()))
    all_rhat = float(az.rhat(phylo_fit).to_array().max())
    min_bulk_ess = float(az.ess(phylo_fit, method="bulk").to_array().min())
    print("All-parameter maximum R-hat:", all_rhat, "| minimum bulk ESS:", min_bulk_ess)
    az.plot_ppc(posterior_check, num_pp_samples=50)
    plt.suptitle(LABEL + " — posterior predictive check")
    plt.show()
    assert all_rhat < 1.01, "Chains have not converged"
    assert min_bulk_ess >= 400, "Insufficient effective samples"
    assert int(phylo_fit.sample_stats.diverging.sum()) == 0, "Resolve divergences"
    draws = phylo_fit.posterior.beta.stack(sample=("chain", "draw"))
    posterior_effects = pd.DataFrame({{
        "coefficient": X_phylo.columns,
        "rate_ratio_median": np.exp(draws.median("sample").values),
        "rate_ratio_q025": np.exp(draws.quantile(.025, dim="sample").values),
        "rate_ratio_q975": np.exp(draws.quantile(.975, dim="sample").values),
    }})
    display(posterior_effects)
    (ROOT / "runs").mkdir(exist_ok=True)
    phylo_fit.to_netcdf(ROOT / "runs" / ("{analysis_name}_phylo_" + MODE + ".nc"))
    posterior_effects.to_csv(ROOT / "runs" / ("{analysis_name}_phylo_effects_" + MODE + ".csv"), index=False)
else:
    print("Posterior MCMC was NOT run. Enable AMP_STUDY_PHYLO_MCMC=1 for posterior estimates.")
''')


common_data_note = '''
The biological unit is a genome. AMP richness counts **distinct peptide sequences
per genome**, after exact parent/source mapping and a complete scoring panel.
The exposure counts the corresponding distinct sequences screened under one
fixed protocol, model version and target panel. Missing or failed screens remain
missing. Repeated parent copies and multiple target scores do not multiply counts.

Use each dataset/cohort separately. ATB's Results report any-target <80 µM, while Methods say median-panel ≤80 µM;
that conflict needs resolution before a real ATB analysis. The demo illustrates
an any-target rule. The archaeome study uses a separate mean-panel ≤100 µM rule. Preserve the
original rule in `bundle.json`; harmonizing thresholds also requires harmonized
screening and source coverage. Predicted MIC and measured MIC are distinct fields.
'''

references = '''
Sources and methods:

- [GenomeSPOT paper](https://doi.org/10.1101/2024.03.22.586313)
- [AllTheBacteria paper](https://doi.org/10.1101/2024.03.08.584059)
- [Archaeome paper](https://doi.org/10.1038/s41564-025-02061-0) and
  [Mendeley data deposit](https://doi.org/10.17632/d8yzgtdrcp.3)
- [AMR annotation documentation](https://allthebacteria.org/docs/amr/)
- [ATB phylogeny documentation](https://allthebacteria.org/docs/phylogeny/)
- [Negative-binomial likelihood](https://www.pymc.io/projects/docs/en/v5.10.4/api/distributions/generated/pymc.NegativeBinomial.html)

See `study/data-contract.md`, `study/analysis-plan.md` and `study/source-search.md`
for input definitions, reproducibility, exclusions and remaining source gaps.
'''

ecology = [
md('''# 1. Ecology and encrypted-AMP richness

Question: do GenomeSPOT temperature, pH, salinity and oxygen-tolerance predictions
explain the rate of high-scoring encrypted AMP sequences after accounting for
screening exposure, proteome size and taxonomy?

**The saved regression outputs use explicitly simulated data.** Actual public
source audits and source-genome predictions are shown separately. Switch to
real mode only with a validated complete screening bundle.
'''), code(SETUP), md("## Available source data"), code(AUDIT), md(common_data_note), code(LOAD),
md('''## Primary model and inclusion criteria

Use quality-screened genomes with successful complete AMP screens and unflagged
predictions for all four primary traits. Keep all exclusion reasons. Fit an NB2
count model with log(screened sequences / 100,000) as a fixed offset. Adjust for
proteome amino-acid count, GC fraction, habitat and family fixed effects. Species
clusters provide robust standard errors; the phylogenetic sensitivity model below
adds branch-length-based covariance. The four ecology coefficients form one
prespecified multiple-testing family (Benjamini–Hochberg correction).

An effect is a rate ratio per cohort standard deviation of the predictor, not a
change in organism growth or antimicrobial activity caused by changing habitat.
'''),
code('''
data, excluded = analysis_subset(cohort, require_traits=True)
traits = [f"trait_{t}" for t in TRAITS]
data, scales = standardized(data, traits + ["gc_fraction", "log_proteome_aa"])
display(pd.Series({"included": len(data), "excluded": len(excluded), "families": data.family.nunique(), "species_clusters": data.species.nunique()}))
display(excluded.head(10))
display(data[traits].corr())
formula = "amp_count ~ " + " + ".join(f"z_{t}" for t in traits) + " + z_gc_fraction + z_log_proteome_aa + C(habitat) + C(family)"
fit = fit_nb(data, formula)
effects = effect_table(fit, [f"z_{t}" for t in traits])
display(effects)
habitat_effects = effect_table(fit, [term for term in fit.params.index if term.startswith("C(habitat)")])
display(habitat_effects)
print("NB2 overdispersion alpha:", fit.params["alpha"])
print("Observed zero fraction:", data.amp_count.eq(0).mean())
print("Predicted zero fraction:", np.mean((1 + fit.params["alpha"] * fit.predict()) ** (-1 / fit.params["alpha"])))
'''),
code('''
fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
axes[0].hist(data.amp_count, bins=25, color="#247b8a")
axes[0].set(xlabel="Distinct qualifying AMP sequences per genome", ylabel="Genomes")
y = np.arange(len(effects))
axes[1].errorbar(effects.rate_ratio, y, xerr=[effects.rate_ratio-effects.ci_low, effects.ci_high-effects.rate_ratio], fmt="o", color="#9c315a")
axes[1].axvline(1, color="grey", linestyle="--")
axes[1].set_yticks(y, effects.term.str.replace("z_trait_", "", regex=False))
axes[1].set(xlabel="Adjusted rate ratio per predictor SD (95% interval)", xscale="log")
from matplotlib.ticker import FuncFormatter, NullLocator
axes[1].set_xticks(np.geomspace(effects.ci_low.min() * .98, effects.ci_high.max() * 1.02, 4))
axes[1].xaxis.set_major_formatter(FuncFormatter(lambda x, position: f"{x:.2g}"))
axes[1].xaxis.set_minor_locator(NullLocator())
fig.suptitle(LABEL)
fig.tight_layout()
save_figure(fig, "ecology_primary")
plt.show()
'''),
md('''## Prediction beyond taxonomy

Compare a numeric size/GC baseline with size/GC plus the four traits, holding
entire families out. Scaling is estimated using each training fold. Family and
habitat fixed effects are omitted from **both prediction models** so unseen
categories do not enter through fitted labels. Lower held-out Poisson deviance
means better mean prediction; it is not a likelihood test or accuracy of the
underlying GenomeSPOT traits. The two models use identical folds.
'''), code('''
cv = grouped_ecology_cv(data)
display(cv)
display(cv.pivot(index="fold", columns="model", values="poisson_deviance"))
'''),
md('''## Sensitivity analyses

Repeat the primary model at a stricter MIC cutoff, and with finite flagged
GenomeSPOT outputs retained. These are sensitivity analyses of the same cohort,
not extra independent discoveries. Missing predictions are never imputed as zero.
If the source table is complete only below the primary threshold, stricter
thresholds are identifiable but looser thresholds are not.
'''), code('''
sensitivity_rows = []
strict_spec = {**spec, "threshold_um": spec["threshold_um"] / 2}
strict_cohort, _, _ = prepare_cohort(tables, strict_spec)
for label, input_cohort, keep_flags in [("stricter_cutoff", strict_cohort, False), ("finite_flagged_included", cohort, True)]:
    subset, _ = analysis_subset(input_cohort, require_traits=not keep_flags)
    subset = subset.dropna(subset=traits).reset_index(drop=True)
    subset, _ = standardized(subset, traits + ["gc_fraction", "log_proteome_aa"])
    result = effect_table(fit_nb(subset, formula), [f"z_{t}" for t in traits])
    result["sensitivity"] = label
    result["n"] = len(subset)
    sensitivity_rows.append(result)
sensitivities = pd.concat(sensitivity_rows, ignore_index=True)
display(sensitivities)
'''),
md('''## Phylogenetic sensitivity

Use the supplied rooted, branch-length tree rather than a taxonomy dendrogram.
The NB2 mixed model includes a Brownian phylogenetic effect and the same numeric
ecology/size/GC covariates plus habitat. Family fixed effects are omitted here
because the tree accounts for ancestry. The covariance, prior predictive model
and likelihood are exercised below. Posterior sampling is opt-in and checked
for R-hat, effective sample size and divergent transitions.
'''), phylo_cell("1 + z_trait_temperature + z_trait_ph + z_trait_salinity + z_trait_oxygen + z_gc_fraction + z_log_proteome_aa + C(habitat)", "ecology"),
export_cell("ecology", {"cohort.csv": "data", "excluded.csv": "excluded", "effects.csv": "effects", "habitat_effects.csv": "habitat_effects", "cross_validation.csv": "cv", "sensitivities.csv": "sensitivities"}),
md('''## Interpretation and limits

Ask whether effects are stable across cutoffs, taxa, flags and family-held-out
prediction. Report effect sizes and intervals, including null results. GenomeSPOT
traits describe predicted growth preferences, not measured habitat conditions.
Protein composition contributes to both physiology and AMP prediction, so a
correlation may reflect shared model inputs. Before a biological claim, add
composition/screen-coverage controls and independent strain-resolved measurements.
Measurement-error sensitivity cannot treat local cross-validation RMSE as a
calibrated per-genome confidence interval.
'''), md(references),
]

evolution = [md('''# 2. Evolution and encrypted-AMP richness

Question: are particular clades enriched in qualifying encrypted AMP sequences
after accounting for screening exposure, proteome size and relatedness?

**Saved inferential results use simulated data.** Real published genus totals
and actual source-genome predictions are descriptive checks with explicit scope.
'''), code(SETUP), code(AUDIT),
md('''## A real published result that must not be overinterpreted

The archaeome paper's Supplementary Table 1 contains the ten genera with the
highest and ten with the lowest reported active-EP counts, along with NCBI gene
totals. This is an outcome-selected genus sample. It cannot support an unbiased
genome-level regression or establish screening effort. Historical names are
preserved; taxonomic harmonization requires a separate explicit mapping.
'''), code('''
published = pd.read_csv(ROOT / "examples/amp-study/published_genus_table.csv")
display(published)
assert len(published) == 20
print("Selected published genus rows only; no significance test is fitted to this table.")
'''), md(common_data_note), code(LOAD),
md('''## Taxonomic coverage and descriptive rates

Count genomes, families and species before looking for enriched clades. Raw
ratios are descriptive; modeling uses counts and exposure offsets. Use an
outcome-independent balanced sampling design and avoid selecting clades for
testing merely because they have the largest observed counts.
'''), code('''
data, excluded = analysis_subset(cohort, require_traits=True)
data, scales = standardized(data, [f"trait_{t}" for t in TRAITS] + ["gc_fraction", "log_proteome_aa"])
clades = data.groupby("phylum").agg(genomes=("genome_id", "size"), species=("species", "nunique"), amp_count=("amp_count", "sum"), screened=("n_screened", "sum"))
clades["amps_per_100000_screened"] = 100000 * clades.amp_count / clades.screened
display(clades)
assert (clades.genomes >= 10).all(), "Prespecify handling of sparsely sampled clades"
ax = clades.amps_per_100000_screened.plot.bar(color="#247b8a", figsize=(7, 3.5))
ax.set(ylabel="Qualifying sequences\\nper 100,000 screened", xlabel="Phylum", title=LABEL + " — descriptive rates")
plt.xticks(rotation=20, ha="right")
plt.tight_layout(); save_figure(plt.gcf(), "evolution_rates"); plt.show()
'''),
md('''## Initial clade model

Compare phyla relative to an explicit reference clade, adjusting for ecology,
proteome size, GC and habitat. Species-clustered intervals are useful as an
initial taxonomic analysis. They do not replace the phylogenetic mixed model.
Adjust the clade contrast family for multiple testing. Repeat at other ranks
only as prespecified sensitivity analyses with appropriate multiplicity control.
'''), code('''
clade_formula = "amp_count ~ C(phylum) + z_log_proteome_aa + z_gc_fraction + z_trait_temperature + z_trait_ph + z_trait_salinity + z_trait_oxygen + C(habitat)"
clade_fit = fit_nb(data, clade_formula)
clade_effects = effect_table(clade_fit, [x for x in clade_fit.params.index if x.startswith("C(phylum)")])
display(clade_effects)
print("Reference phylum:", sorted(data.phylum.unique())[0])
'''), md('''## Primary relatedness model

Let `C` be Brownian shared-branch covariance from the explicitly rooted tree,
normalized to a unit diagonal. Fit `log(mu) = log(exposure/100000) + X beta + u`,
where `u ~ Normal(0, sigma_phylo² C)` and counts follow a negative binomial.
`X` includes the clade contrasts, proteome size, GC, ecology and habitat. Clade
effects and the phylogenetic effect can be weakly identified; report uncertainty
and prior sensitivity. A taxonomy label, cluster-robust standard error or
dendrogram alone is not a phylogenetic adjustment.

Run bacteria and archaea separately with appropriate marker trees. The ATB
bacterial marker tree is not a valid archaeal tree. This implementation checks
unique tips, exact sample alignment and finite nonnegative branch lengths, and
limits dense covariance to 500 genomes. Root choice is an assumption to test.
'''), phylo_cell("1 + C(phylum) + z_log_proteome_aa + z_gc_fraction + z_trait_temperature + z_trait_ph + z_trait_salinity + z_trait_oxygen + C(habitat)", "evolution"),
code('''
fig, ax = plt.subplots(figsize=(5, 4))
image = ax.imshow(phylo_cov[:80, :80], vmin=0, vmax=1, cmap="viridis")
fig.colorbar(image, ax=ax, label="Shared-branch correlation")
ax.set(title=LABEL + " — phylogenetic covariance", xlabel="Aligned genome index (first 80)", ylabel="Aligned genome index (first 80)")
plt.tight_layout(); save_figure(fig, "evolution_covariance"); plt.show()
'''),
export_cell("evolution", {"cohort.csv": "data", "excluded.csv": "excluded", "clade_coverage.csv": "clades.reset_index()", "clade_effects.csv": "clade_effects"}),
md('''## Evidence needed for an evolutionary conclusion

Replicate the clade comparison using an outcome-independent genome panel and
alternative plausible trees/rooting. Compare rates at equal screened sequence
coverage. Repeat within clades to check whether one species or heavily sequenced
lineage drives the result. Use posterior predictive checks and sampling
diagnostics before reporting phylogenetic posterior estimates. A clade association
alone does not demonstrate selection, ecological function or a particular gene
transfer mechanism. Published curated-protein coverage can vary strongly by clade.
'''), md(references)]

resistance = [md('''# 3. Resistance burden and encrypted-AMP opportunity

Question: do ATB genomes with more annotated AMR determinants have greater AMP
richness, stronger predicted activity among qualifying candidates, or greater
candidate-cluster diversity at equal sampling effort?

**Saved regression outputs use simulated data.** This notebook tests a framework,
not an observed AMR–AMP or co-evolution relationship. Genome AMR burden refers to
the source organism, whereas an AMP's target panel refers to organisms against
which the predictor estimates activity. These are different biological roles.
'''), code(SETUP), code(AUDIT), md(common_data_note), code(LOAD),
md('''## Build resistance burden without false zeros

Use AMRFinderPlus element type `AMR` only; exclude stress/virulence annotations.
Count distinct gene symbols as the primary burden and distinct reported drug
classes as a sensitivity measure. Retain locus IDs to identify duplicate calls.
A successful AMR run without hits contributes zero. Failed, not-run and absent
status entries contribute missing values. Do not interpret these genetic calls
as measured susceptibility, acquired resistance only, or interchangeable effects
of every determinant. An acquired-only analysis needs an independently curated
intrinsic/acquired classification and a fixed database release.
'''), code('''
if COHORT != "atb":
    raise ValueError("This AMR analysis requires an ATB cohort with matched AMRFinderPlus results")
data, excluded = analysis_subset(cohort, require_traits=True, require_amr=True)
columns = [f"trait_{t}" for t in TRAITS] + ["gc_fraction", "log_proteome_aa", "amr_gene_count", "amr_class_count"]
data, scales = standardized(data, columns)
display(data[["genome_id", "species", "amr_gene_count", "amr_class_count", "amp_count", "n_screened"]].head())
display(cohort.amr_status.value_counts(dropna=False))
'''),
md('''## Richness across and within species

Model counts with the screening offset, proteome size, GC, ecology and habitat.
The across-cohort model includes family effects. The within-species model uses
species fixed effects and requires actual AMR variation within at least ten
species. This helps distinguish a lineage's baseline biology from isolate-level
co-occurrence. It still cannot establish co-evolution or causality.
'''), code('''
covariates = "z_log_proteome_aa + z_gc_fraction + z_trait_temperature + z_trait_ph + z_trait_salinity + z_trait_oxygen + C(habitat)"
across_formula = "amp_count ~ z_amr_gene_count + " + covariates + " + C(family)"
across_fit = fit_nb(data, across_formula)
across_effect = effect_table(across_fit, ["z_amr_gene_count"])
across_effect["comparison"] = "across cohort; family adjusted"
variable_species = data.groupby("species").amr_gene_count.nunique()
variable_species = variable_species[variable_species >= 2].index
within = data[data.species.isin(variable_species)].copy()
assert within.species.nunique() >= 10, "Too few species with within-species AMR variation"
within_formula = "amp_count ~ z_amr_gene_count + " + covariates + " + C(species)"
within_fit = fit_nb(within, within_formula)
within_effect = effect_table(within_fit, ["z_amr_gene_count"])
within_effect["comparison"] = "within species"
richness_effects = pd.concat([across_effect, within_effect], ignore_index=True)
display(richness_effects)
# Alternative burden definition is a prespecified sensitivity, not a second claim.
class_fit = fit_nb(data, across_formula.replace("z_amr_gene_count", "z_amr_class_count"))
display(effect_table(class_fit, ["z_amr_class_count"]))
'''),
md('''## Conditional potency and equal-effort diversity

Use each genome's median log10 predicted MIC among qualifying unique candidates,
requiring at least five candidates. Lower values imply stronger **predicted**
activity under the declared any-target or mean-panel rule. This is conditional
on passing the threshold; a hit-only catalog cannot estimate the potency of all
fragments. Include log(candidate count) to reduce count-dependent selection bias.
Never substitute a best-of-many MIC as the main potency endpoint.

For distinctness, use cohort-wide, independently generated peptide-cluster IDs
and analytical rarefaction to ten candidate sequences per genome. This controls
sampling count, but it does not prove novel sequence mechanisms. Missing clusters
leave diversity unavailable. Experimental MICs need a separate analysis preserving
target identities, assay conditions and left/right censoring; this code does not
replace censored values with an arbitrary numeric MIC.
'''), code('''
import statsmodels.formula.api as smf
potency = data[data.n_potency_candidates.ge(5)].copy()
potency["log_amp_count"] = np.log(potency.amp_count)
potency = potency[potency.species.isin(variable_species)].copy()
assert len(potency) >= 100 and potency.species.nunique() >= 10
potency_formula = "median_predicted_log10_mic_um ~ z_amr_gene_count + log_amp_count + " + covariates + " + C(species)"
potency_fit = smf.ols(potency_formula, potency).fit(cov_type="cluster", cov_kwds={"groups": potency.species})
potency_effect = pd.DataFrame({"coefficient": [potency_fit.params["z_amr_gene_count"]], "ci_low": [potency_fit.conf_int().loc["z_amr_gene_count", 0]], "ci_high": [potency_fit.conf_int().loc["z_amr_gene_count", 1]], "p_value": [potency_fit.pvalues["z_amr_gene_count"]], "n_genomes": [len(potency)]})
display(potency_effect)
print("Multiplicative change in conditional predicted MIC per AMR SD:", 10 ** potency_effect.coefficient.iloc[0])
diversity = data.dropna(subset=["rarefied_cluster_richness_10"]).copy()
diversity["log_amp_count"] = np.log(diversity.amp_count)
diversity_model = None
if len(diversity) >= 100 and diversity.species.nunique() >= 10:
    diversity_formula = "rarefied_cluster_richness_10 ~ z_amr_gene_count + log_amp_count + " + covariates + " + C(family)"
    diversity_model = smf.ols(diversity_formula, diversity).fit(cov_type="cluster", cov_kwds={"groups": diversity.species})
    display(pd.DataFrame({"amr_coefficient": [diversity_model.params["z_amr_gene_count"]], "p_value": [diversity_model.pvalues["z_amr_gene_count"]], "n_genomes": [len(diversity)]}))
else:
    print("Cluster diversity: insufficient mapped clusters/qualifying genomes; result unavailable.")
# Three distinct endpoints form a testing family; report joint FDR adjustment.
from statsmodels.stats.multitest import multipletests
endpoint_p = {"within_species_richness": within_fit.pvalues["z_amr_gene_count"], "conditional_potency": potency_fit.pvalues["z_amr_gene_count"], "rarefied_diversity": diversity_model.pvalues["z_amr_gene_count"] if diversity_model is not None else np.nan}
endpoint_tests = pd.DataFrame({"endpoint": list(endpoint_p), "p_value": list(endpoint_p.values())})
endpoint_tests["q_value_bh"] = multipletests(endpoint_tests.p_value.fillna(1), method="fdr_bh")[1]
endpoint_tests.loc[endpoint_tests.p_value.isna(), "q_value_bh"] = np.nan
display(endpoint_tests)
'''),
code('''
fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
axes[0].scatter(data.amr_gene_count, 100000 * data.amp_count / data.n_screened, alpha=.45, s=16, color="#247b8a")
axes[0].set(xlabel="Distinct annotated AMR genes", ylabel="Qualifying AMP sequences\\nper 100,000 screened")
axes[1].scatter(potency.amr_gene_count, potency.median_predicted_log10_mic_um, alpha=.5, s=16, color="#9c315a")
axes[1].set(xlabel="Distinct annotated AMR genes", ylabel="Conditional median log10\\npredicted MIC (µM)")
fig.suptitle(LABEL + " — descriptive plots; adjusted models above")
fig.tight_layout(); save_figure(fig, "resistance_descriptive"); plt.show()
'''),
md("## Phylogenetic sensitivity of the richness association\n\nUse the same rooted-tree NB2 model with AMR burden as a predictor. Report disagreement between across-species, within-species and phylogenetic results rather than choosing the largest effect."),
phylo_cell("1 + z_amr_gene_count + z_log_proteome_aa + z_gc_fraction + z_trait_temperature + z_trait_ph + z_trait_salinity + z_trait_oxygen + C(habitat)", "resistance"),
export_cell("resistance", {"cohort.csv": "data", "excluded.csv": "excluded", "richness_effects.csv": "richness_effects", "conditional_potency.csv": "potency_effect", "endpoint_tests.csv": "endpoint_tests"}),
md('''## Interpretation and follow-up

Co-occurrence can arise from taxonomy, mobile elements, clinical sampling,
annotation coverage and common sequence composition. A co-evolution argument
needs additional evidence such as concordant gains/losses on a suitable tree,
gene-context or transfer analysis and independent validation. Do not infer
natural expression, peptide release, ecological antagonism or clinical utility
from these predicted scores. Repeat with curated acquired-only AMR calls,
alternative database releases and measured susceptibility where available.
'''), md(references)]


for filename, cells in [("01_ecology_amp_richness.ipynb", ecology), ("02_evolution_amp_richness.ipynb", evolution), ("03_resistance_amp_opportunity.ipynb", resistance)]:
    notebook = nbf.v4.new_notebook(cells=cells, metadata={"kernelspec": {"display_name": "Python 3 (GenomeSPOT study)", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.11"}})
    nbf.validate(notebook)
    nbf.write(notebook, ROOT / "notebooks" / filename)
    print(filename, len(cells), "cells")
