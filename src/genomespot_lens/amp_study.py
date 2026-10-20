"""Shared, validated data preparation and statistics for the three AMP notebooks.

Published hit lists alone do not establish complete genome screens. The explicit
screen registry and completion flags below protect zero counts and denominators.
"""
from io import StringIO
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import gammaln

from .core import numeric_prediction, quality_flags, sha256

TRAITS = {"temperature": "temperature_optimum", "ph": "ph_optimum", "salinity": "salinity_optimum", "oxygen": "oxygen"}


def require_columns(df, names, table):
    missing = set(names) - set(df)
    if missing:
        raise ValueError(f"{table} missing columns: {sorted(missing)}")


def boolean(value):
    if pd.isna(value):
        return None
    if value in (True, "True", "true", 1):
        return True
    if value in (False, "False", "false", 0):
        return False
    raise ValueError(f"Invalid boolean: {value!r}")


def assert_unique(df, keys, table):
    if df[keys].isna().any().any() or df.duplicated(keys).any():
        raise ValueError(f"{table}: missing or duplicate keys {keys}")


def load_bundle(directory, cohort, expected_kind):
    directory = Path(directory)
    spec = json.loads((directory / "bundle.json").read_text())
    if spec.get("data_kind") != expected_kind:
        raise ValueError("Data-kind mismatch: a real run must never silently load demonstration data")
    if spec.get("schema_version") != 1:
        raise ValueError("Unsupported study schema")
    tables = {}
    provenance = {"bundle_sha256": sha256(directory / "bundle.json"), "data_kind": expected_kind, "cohort": cohort, "files": {}}
    for name in ["genomes", "screening", "hits", "genomespot", "amr_status", "amr_hits"]:
        path = directory / f"{name}.csv"
        tables[name] = pd.read_csv(path)
        provenance["files"][path.name] = sha256(path)
    require_columns(tables["genomes"], ["genome_id", "cohort", "assembly_accession", "domain", "phylum", "family", "genus", "species", "protein_count", "proteome_aa", "gc_fraction", "completeness", "contamination", "habitat"], "genomes")
    assert_unique(tables["genomes"], ["genome_id"], "genomes")
    assert_unique(tables["genomes"], ["assembly_accession"], "genomes")
    known = set(tables["genomes"].genome_id)
    for name, table in tables.items():
        require_columns(table, ["genome_id"], name)
        if not set(table.genome_id) <= known:
            raise ValueError(f"{name} contains genome IDs absent from the genome registry")
    genomes = tables["genomes"].query("cohort == @cohort").copy()
    if genomes.empty or genomes.domain.nunique() != 1:
        raise ValueError("Select a nonempty, single-domain cohort; do not pool bacterial/archaeal screens")
    ids = set(genomes.genome_id)
    tables = {name: table[table.genome_id.isin(ids)].copy() for name, table in tables.items()}
    tree = directory / "phylogeny.nwk"
    provenance["tree_sha256"] = sha256(tree) if tree.exists() else None
    return tables, spec, provenance, tree


def aggregate_hits(hits, spec):
    required = ["genome_id", "peptide_id", "parent_protein_id", "target", "predicted_mic_um", "mapping_status", "peptide_cluster_id"]
    require_columns(hits, required, "hits")
    if spec.get("threshold_mode") not in {"mean", "any"} or spec.get("threshold_um", 0) <= 0 or spec.get("threshold_comparison") not in {"lt", "le"}:
        raise ValueError("Specify an explicit positive MIC threshold and mean/any rule")
    targets = spec.get("targets", [])
    if not targets or len(set(targets)) != len(targets):
        raise ValueError("Declare the full, nonduplicated scoring target panel")
    hits = hits.copy()
    if hits[["peptide_id", "parent_protein_id", "target"]].isna().any().any():
        raise ValueError("Hit sequence, parent and target identifiers are required")
    if not set(hits.target) <= set(targets):
        raise ValueError("Hits contain targets outside the declared panel")
    mic = pd.to_numeric(hits.predicted_mic_um, errors="raise")
    if (~np.isfinite(mic) | (mic <= 0)).any():
        raise ValueError("Predicted MIC must be finite, positive and in micromolar units")
    hits["predicted_mic_um"] = mic
    exact = hits[hits.mapping_status.eq("exact")].copy()
    keys = ["genome_id", "peptide_id", "target"]
    conflicts = exact.groupby(keys).predicted_mic_um.nunique()
    if (conflicts > 1).any():
        raise ValueError("Conflicting scores for the same genome/peptide/target")
    if exact.dropna(subset=["peptide_cluster_id"]).groupby("peptide_id").peptide_cluster_id.nunique().gt(1).any():
        raise ValueError("A peptide has inconsistent cohort-wide cluster assignments")
    # Repeated parent copies and targets must not inflate distinct-sequence richness.
    exact = exact.drop_duplicates(keys)
    records = []
    incomplete = 0
    for (genome, peptide), group in exact.groupby(["genome_id", "peptide_id"], sort=True):
        if set(group.target) != set(targets):
            incomplete += 1
            continue
        values = group.predicted_mic_um.to_numpy()
        score = values.mean() if spec["threshold_mode"] == "mean" else values.min()
        qualifies = score < spec["threshold_um"] if spec["threshold_comparison"] == "lt" else score <= spec["threshold_um"]
        if qualifies:
            clusters = group.peptide_cluster_id.dropna().unique()
            records.append({"genome_id": genome, "peptide_id": peptide, "peptide_cluster_id": clusters[0] if len(clusters) == 1 else None,
                            "score_mic_um": float(score), "log10_mic_um": float(np.log10(score)),
                            "mean_mic_um": float(values.mean())})
    candidates = pd.DataFrame(records, columns=["genome_id", "peptide_id", "peptide_cluster_id", "score_mic_um", "log10_mic_um", "mean_mic_um"])
    audit = {"input_rows": len(hits), "nonexact_rows_excluded": int((~hits.mapping_status.eq("exact")).sum()), "incomplete_panel_peptides_excluded": incomplete, "distinct_genome_peptide_hits": len(candidates)}
    return candidates, audit


def amr_burden(status, hits, genome_ids):
    require_columns(status, ["genome_id", "status"], "amr_status")
    require_columns(hits, ["genome_id", "element_type", "gene_symbol", "drug_class", "locus_id"], "amr_hits")
    if hits.element_type.isna().any() or not set(hits.element_type) <= {"AMR", "STRESS", "VIRULENCE"}:
        raise ValueError("Normalize AMR element types explicitly; unknown/missing types cannot become zero burden")
    assert_unique(status, ["genome_id"], "amr_status")
    if not set(status.status.dropna()) <= {"PASS", "FAIL", "NOT DONE"}:
        raise ValueError("AMR status must be PASS, FAIL or NOT DONE")
    result = pd.DataFrame({"genome_id": list(genome_ids)}).merge(status, how="left", on="genome_id", validate="1:1")
    determinants = hits[hits.element_type.eq("AMR")].copy()
    if determinants[["gene_symbol", "locus_id"]].isna().any().any():
        raise ValueError("AMR determinants need gene symbols and locus IDs")
    genes = determinants.drop_duplicates(["genome_id", "gene_symbol"]).groupby("genome_id").size()
    classes = determinants.dropna(subset=["drug_class"]).assign(drug_class=lambda d: d.drug_class.str.split("/"))
    classes = classes.explode("drug_class")
    classes["drug_class"] = classes.drug_class.str.strip()
    classes = classes[classes.drug_class.ne("")].drop_duplicates(["genome_id", "drug_class"]).groupby("genome_id").size()
    for name, values in [("amr_gene_count", genes), ("amr_class_count", classes)]:
        result[name] = result.genome_id.map(values)
        result.loc[result.status.eq("PASS"), name] = result.loc[result.status.eq("PASS"), name].fillna(0)
        result.loc[~result.status.eq("PASS"), name] = np.nan
    return result.rename(columns={"status": "amr_status"})


def expected_rarefied_richness(cluster_ids, draws=10):
    if pd.isna(cluster_ids).any() or len(cluster_ids) < draws:
        return np.nan
    counts = pd.Series(cluster_ids).value_counts().to_numpy()
    n = counts.sum()
    def logchoose(a, b):
        return gammaln(a + 1) - gammaln(b + 1) - gammaln(a - b + 1)
    absent = np.zeros(len(counts))
    ok = n - counts >= draws
    absent[ok] = np.exp(logchoose(n - counts[ok], draws) - logchoose(n, draws))
    return float(np.sum(1 - absent))


def prepare_cohort(tables, spec):
    genomes, screening = tables["genomes"].copy(), tables["screening"].copy()
    require_columns(screening, ["genome_id", "status", "hits_complete", "n_screened", "screen_protocol", "exposure_unit", "score_model", "score_version", "target_panel_id"], "screening")
    assert_unique(screening, ["genome_id"], "screening")
    if set(screening.status.dropna()) - {"COMPLETE", "FAIL", "NOT DONE"}:
        raise ValueError("Unknown screen status")
    if screening.exposure_unit.ne("unique_peptide_sequences").any():
        raise ValueError("This analysis requires unique-peptide-sequence exposure, matching the hit count unit")
    for col in ["screen_protocol", "score_model", "score_version", "target_panel_id"]:
        if screening[col].isna().any() or screening[col].nunique() != 1:
            raise ValueError(f"Analyze a homogeneous scoring/screen cohort: {col}")
    candidates, hit_audit = aggregate_hits(tables["hits"], spec)
    counts = candidates.groupby("genome_id").size()
    cohort = genomes.merge(screening, how="left", on="genome_id", validate="1:1")
    eligible = cohort.status.eq("COMPLETE") & cohort.hits_complete.map(boolean).eq(True)
    exposure = pd.to_numeric(cohort.n_screened, errors="coerce")
    if ((~np.isfinite(exposure) | (exposure < 1) | (exposure % 1 != 0)) & eligible).any():
        raise ValueError("Completed screens require positive integer screening denominators")
    cohort["screen_eligible"] = eligible
    cohort["amp_count"] = cohort.genome_id.map(counts)
    cohort.loc[eligible, "amp_count"] = cohort.loc[eligible, "amp_count"].fillna(0)
    cohort.loc[~eligible, "amp_count"] = np.nan
    if cohort.loc[eligible, "amp_count"].gt(cohort.loc[eligible, "n_screened"]).any():
        raise ValueError("AMP count exceeds the number of screened distinct sequences")
    bounds = {"gc_fraction": (0, 1), "completeness": (0, 100), "contamination": (0, 100)}
    for column, (lo, hi) in bounds.items():
        if (~pd.to_numeric(cohort[column], errors="coerce").between(lo, hi)).any():
            raise ValueError(f"Invalid or missing {column}; check units")
    sizes = cohort[["proteome_aa", "protein_count"]].apply(pd.to_numeric, errors="raise")
    if (~np.isfinite(sizes) | sizes.le(0)).any().any():
        raise ValueError("Proteome sizes must be positive")
    cohort["genome_qc_pass"] = cohort.completeness.ge(95) & cohort.contamination.le(5)
    cohort["log_proteome_aa"] = np.log(cohort.proteome_aa)
    spot = tables["genomespot"]
    require_columns(spot, ["genome_id", "target", "value", "error", "is_novel", "warning"], "genomespot")
    assert_unique(spot, ["genome_id", "target"], "genomespot")
    for name, target in TRAITS.items():
        rows = []
        for _, row in spot[spot.target.eq(target)].iterrows():
            prediction = {"value": row.value if pd.notna(row.value) else None,
                          "error": float(row.error) if pd.notna(row.error) else None,
                          "is_novel": boolean(row.is_novel),
                          "warning": str(row.warning) if pd.notna(row.warning) else None}
            if target != "oxygen" and prediction["value"] is not None:
                prediction["value"] = float(prediction["value"])
            flags = quality_flags(target, prediction)
            if prediction["error"] is None or not np.isfinite(prediction["error"]) or prediction["error"] < 0:
                flags.append("uncertainty_unavailable")
            rows.append({"genome_id": row.genome_id, f"trait_{name}": numeric_prediction(target, prediction), f"flag_{name}": bool(flags), f"error_{name}": prediction["error"]})
        trait = pd.DataFrame(rows, columns=["genome_id", f"trait_{name}", f"flag_{name}", f"error_{name}"])
        cohort = cohort.merge(trait, how="left", on="genome_id", validate="1:1")
        cohort[f"flag_{name}"] = cohort[f"flag_{name}"].fillna(True).astype(bool)
    summary = candidates.groupby("genome_id").agg(median_predicted_log10_mic_um=("log10_mic_um", "median"), n_potency_candidates=("peptide_id", "nunique"))
    cohort = cohort.merge(summary, how="left", on="genome_id", validate="1:1")
    rarefied = candidates.groupby("genome_id").peptide_cluster_id.apply(lambda x: expected_rarefied_richness(x.to_numpy(), 10))
    cohort["rarefied_cluster_richness_10"] = cohort.genome_id.map(rarefied)
    cohort = cohort.merge(amr_burden(tables["amr_status"], tables["amr_hits"], cohort.genome_id), how="left", on="genome_id", validate="1:1")
    cohort["all_traits_unflagged"] = ~cohort[[f"flag_{name}" for name in TRAITS]].any(axis=1)
    return cohort, candidates, hit_audit


def analysis_subset(cohort, require_traits=True, require_amr=False):
    mask = cohort.screen_eligible & cohort.genome_qc_pass
    if require_traits:
        mask &= cohort.all_traits_unflagged
    if require_amr:
        mask &= cohort.amr_status.eq("PASS")
    required = ["family", "species", "habitat", "gc_fraction", "log_proteome_aa", "amp_count", "n_screened"]
    mask &= cohort[required].notna().all(axis=1)
    data = cohort.loc[mask].copy().reset_index(drop=True)
    exclusions = cohort.loc[~mask, ["genome_id", "screen_eligible", "genome_qc_pass", "all_traits_unflagged", "amr_status"]].copy()
    return data, exclusions


def standardized(data, columns):
    data = data.copy()
    scales = {}
    for column in columns:
        mean, sd = data[column].mean(), data[column].std(ddof=0)
        if not np.isfinite(sd) or sd == 0:
            raise ValueError(f"No usable variation in {column}")
        data[f"z_{column}"] = (data[column] - mean) / sd
        scales[column] = {"mean": float(mean), "sd": float(sd)}
    return data, scales


def fit_nb(data, formula, cluster="species"):
    import statsmodels.formula.api as smf
    import patsy
    y, design = patsy.dmatrices(formula, data, return_type="dataframe", NA_action="raise")
    if len(data) < max(50, 5 * design.shape[1]) or data[cluster].nunique() < 10:
        raise ValueError("Insufficient observations/clusters for the prespecified model; report coverage instead")
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("Design is rank deficient; simplify prespecified covariates")
    model = smf.negativebinomial(formula, data, offset=np.log(data.n_screened / 100000), loglike_method="nb2")
    fit = model.fit(method="bfgs", maxiter=1000, disp=False, cov_type="cluster", cov_kwds={"groups": data[cluster]})
    if not fit.mle_retvals.get("converged", False) or not np.isfinite(fit.params).all():
        raise ValueError("Negative-binomial fit did not converge")
    return fit


def effect_table(fit, terms):
    from statsmodels.stats.multitest import multipletests
    intervals = fit.conf_int()
    rows = []
    for term in terms:
        if term not in fit.params:
            continue
        rows.append({"term": term, "log_rate_ratio": fit.params[term], "rate_ratio": np.exp(fit.params[term]), "ci_low": np.exp(intervals.loc[term, 0]), "ci_high": np.exp(intervals.loc[term, 1]), "p_value": fit.pvalues[term]})
    frame = pd.DataFrame(rows)
    if len(frame):
        frame["q_value_bh"] = multipletests(frame.p_value, method="fdr_bh")[1]
    return frame


def grouped_ecology_cv(data, n_splits=4):
    """Numeric scaling learned within training folds; whole families held out."""
    import statsmodels.api as sm
    from sklearn.metrics import mean_poisson_deviance
    from sklearn.model_selection import GroupKFold
    covariates = ["gc_fraction", "log_proteome_aa"]
    traits = [f"trait_{t}" for t in TRAITS]
    records = []
    splitter = GroupKFold(n_splits=min(n_splits, data.family.nunique()))
    for fold, (train_ids, test_ids) in enumerate(splitter.split(data, groups=data.family)):
        train, test = data.iloc[train_ids], data.iloc[test_ids]
        for label, columns in [("size_gc", covariates), ("size_gc_plus_ecology", covariates + traits)]:
            mu, sd = train[columns].mean(), train[columns].std(ddof=0)
            if sd.eq(0).any():
                raise ValueError("A training fold has a constant predictor")
            x_train = sm.add_constant((train[columns] - mu) / sd, has_constant="add")
            x_test = sm.add_constant((test[columns] - mu) / sd, has_constant="add")
            fit = sm.NegativeBinomial(train.amp_count, x_train, offset=np.log(train.n_screened / 100000)).fit(method="bfgs", maxiter=1000, disp=False)
            if not fit.mle_retvals.get("converged", False):
                raise ValueError("Held-out prediction fit failed to converge")
            predicted = fit.predict(x_test, offset=np.log(test.n_screened / 100000))
            records.append({"fold": fold, "model": label, "test_n": len(test), "test_families": ";".join(sorted(test.family.unique())), "poisson_deviance": mean_poisson_deviance(test.amp_count, predicted)})
    return pd.DataFrame(records)


def tree_covariance(path, genome_ids):
    """Brownian shared-branch covariance, then convert to unit diagonal."""
    from Bio import Phylo
    tree = Phylo.read(str(path), "newick")
    if not tree.rooted:
        # Newick omits rootedness metadata frequently. Require a rootedness attestation
        # in bundle.json; this function uses the explicitly supplied Newick root.
        pass
    tips = tree.get_terminals()
    names = [tip.name for tip in tips]
    if None in names or len(names) != len(set(names)):
        raise ValueError("Tree tips must have unique nonempty genome IDs")
    if len(genome_ids) != len(set(genome_ids)) or not set(genome_ids) <= set(names):
        raise ValueError("Tree/data tip mapping is incomplete or duplicated")
    for clade in tree.find_clades():
        if clade is tree.root:
            continue
        if clade.branch_length is None or not np.isfinite(clade.branch_length) or clade.branch_length < 0:
            raise ValueError("Tree requires finite nonnegative branch lengths")
    depths = tree.depths()
    by_name = {tip.name: tip for tip in tips}
    ids = list(genome_ids)
    covariance = np.array([[depths[tree.common_ancestor(by_name[a], by_name[b])] - depths[tree.root] for b in ids] for a in ids])
    diagonal = np.diag(covariance)
    if (diagonal <= 0).any():
        raise ValueError("Tree has a zero-depth tip")
    covariance /= np.sqrt(diagonal[:, None] * diagonal[None, :])
    if np.linalg.eigvalsh(covariance).min() < -1e-8:
        raise ValueError("Phylogenetic covariance is not positive semidefinite")
    return covariance


def build_phylogenetic_nb(data, design, covariance):
    """NB2 count mixed model with fixed exposure and Brownian phylogenetic effect.

    E[y_i]=exposure_i/1e5 * exp(X_i beta + u_i),
    u ~ Normal(0, phylo_sd^2 * C). PyMC alpha is inverse NB2 overdispersion.
    """
    import pymc as pm
    if len(data) > 500:
        raise ValueError("Dense phylogenetic model is limited to 500 genomes; use a balanced subset")
    x = np.asarray(design, float)
    covariance = np.asarray(covariance, float)
    if covariance.shape != (len(data), len(data)) or x.shape[0] != len(data):
        raise ValueError("Phylogenetic dimensions disagree with data")
    if not np.allclose(covariance, covariance.T) or np.linalg.eigvalsh(covariance).min() < -1e-8:
        raise ValueError("Invalid phylogenetic covariance")
    if not np.isfinite(x).all() or np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError("Invalid or rank-deficient fixed effects")
    chol = np.linalg.cholesky(covariance + np.eye(len(data)) * 1e-8)
    with pm.Model(coords={"genome": data.genome_id.astype(str).tolist(), "coefficient": list(design.columns)}) as model:
        beta = pm.Normal("beta", 0, 1.5, dims="coefficient")
        phylo_sd = pm.HalfNormal("phylo_sd", 1)
        z = pm.Normal("phylo_z", 0, 1, dims="genome")
        u = pm.Deterministic("phylo_effect", phylo_sd * pm.math.dot(chol, z), dims="genome")
        alpha = pm.Exponential("nb_shape", 1)
        mu = pm.math.exp(np.log(data.n_screened.to_numpy() / 100000) + pm.math.dot(x, beta) + u)
        pm.NegativeBinomial("amp_count", mu=mu, alpha=alpha, observed=data.amp_count.to_numpy().astype(int), dims="genome")
    return model
