import copy
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("statsmodels")

from genomespot_lens.amp_study import (
    TRAITS, aggregate_hits, amr_burden, analysis_subset, build_phylogenetic_nb,
    expected_rarefied_richness, fit_nb, load_bundle, prepare_cohort,
    standardized, tree_covariance,
)

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def bundle():
    return load_bundle(ROOT / "examples/amp-study/demo", "atb", "demo")


def test_demo_never_silently_used_as_real():
    with pytest.raises(ValueError, match="Data-kind mismatch"):
        load_bundle(ROOT / "examples/amp-study/demo", "atb", "real")


def test_failed_screens_and_amr_remain_missing(bundle):
    tables, spec, _, _ = bundle
    cohort, _, _ = prepare_cohort(tables, spec)
    assert cohort.amp_count.eq(0).any()
    assert cohort.loc[~cohort.screen_eligible, "amp_count"].isna().all()
    assert cohort.loc[~cohort.amr_status.eq("PASS"), "amr_gene_count"].isna().all()
    assert cohort.loc[cohort.amr_status.eq("PASS"), "amr_gene_count"].notna().all()


def test_duplicate_targets_and_parent_copies_do_not_inflate_richness(bundle):
    tables, spec, _, _ = bundle
    expected, _ = aggregate_hits(tables["hits"], spec)
    repeated = tables["hits"].iloc[:10].copy()
    repeated["parent_protein_id"] = "different_parent_copy"
    actual, _ = aggregate_hits(pd.concat([tables["hits"], repeated]), spec)
    pd.testing.assert_frame_equal(actual, expected)


def test_ambiguous_mapping_and_incomplete_panel_are_excluded(bundle):
    tables, spec, _, _ = bundle
    hits = tables["hits"].iloc[:4].copy()
    hits.loc[hits.index[0], "mapping_status"] = "ambiguous"
    result, audit = aggregate_hits(hits, spec)
    assert len(result) == 1
    assert audit["incomplete_panel_peptides_excluded"] == 1
    assert audit["nonexact_rows_excluded"] == 1


def test_conflicting_prediction_is_an_error(bundle):
    tables, spec, _, _ = bundle
    row = tables["hits"].iloc[:1].copy()
    row["predicted_mic_um"] *= 2
    with pytest.raises(ValueError, match="Conflicting"):
        aggregate_hits(pd.concat([tables["hits"], row]), spec)


def test_original_threshold_comparison_is_preserved(bundle):
    tables, spec, _, _ = bundle
    hits = tables["hits"].iloc[:2].copy()
    hits["predicted_mic_um"] = 80
    assert aggregate_hits(hits, {**spec, "threshold_comparison": "lt"})[0].empty
    assert len(aggregate_hits(hits, {**spec, "threshold_comparison": "le"})[0]) == 1


def test_unknown_exposure_or_mixed_score_protocol_rejected(bundle):
    tables, spec, _, _ = bundle
    for col, value in [("exposure_unit", "protein_count"), ("score_version", "different")]:
        altered = {k: v.copy() for k, v in tables.items()}
        altered["screening"].loc[1, col] = value
        with pytest.raises(ValueError):
            prepare_cohort(altered, spec)


def test_amr_only_successful_empty_runs_are_zero():
    status = pd.DataFrame({"genome_id": ["a", "b", "c"], "status": ["PASS", "PASS", "FAIL"]})
    hits = pd.DataFrame([
        {"genome_id": "a", "element_type": "AMR", "gene_symbol": "gene1", "drug_class": "classA/classB", "locus_id": "locus1"},
        {"genome_id": "a", "element_type": "AMR", "gene_symbol": "gene1", "drug_class": "classA/classB", "locus_id": "locus2"},
        {"genome_id": "a", "element_type": "STRESS", "gene_symbol": "other", "drug_class": "stress", "locus_id": "locus3"},
    ])
    actual = amr_burden(status, hits, ["a", "b", "c", "d"]).set_index("genome_id")
    assert actual.loc["a", "amr_gene_count"] == 1
    assert actual.loc["a", "amr_class_count"] == 2
    assert actual.loc["b", "amr_gene_count"] == 0
    assert actual.loc[["c", "d"], "amr_gene_count"].isna().all()


def test_rarefaction_matches_exact_small_example():
    # Of three equally likely pairs from A,A,B, two contain both clusters.
    assert expected_rarefied_richness(np.array(["a", "a", "b"]), 2) == pytest.approx(5 / 3)
    assert np.isnan(expected_rarefied_richness(np.array(["a"]), 2))
    assert np.isnan(expected_rarefied_richness(np.array(["a", None]), 2))


def test_phylogenetic_tip_order_and_shared_branches(tmp_path):
    path = tmp_path / "tree.nwk"
    path.write_text("((a:1,b:1):1,c:2):9;")
    covariance = tree_covariance(path, ["c", "b", "a"])
    np.testing.assert_allclose(covariance, [[1, 0, 0], [0, 1, .5], [0, .5, 1]])
    with pytest.raises(ValueError, match="incomplete"):
        tree_covariance(path, ["unknown"])


def test_branch_lengths_are_required(tmp_path):
    path = tmp_path / "tree.nwk"
    path.write_text("(a,b);")
    with pytest.raises(ValueError, match="branch lengths"):
        tree_covariance(path, ["a", "b"])


def test_screen_exposure_changes_intercept_not_trait_effect(bundle):
    tables, spec, _, _ = bundle
    cohort, _, _ = prepare_cohort(tables, spec)
    data, _ = analysis_subset(cohort)
    data, _ = standardized(data, ["trait_temperature", "gc_fraction", "log_proteome_aa"])
    formula = "amp_count ~ z_trait_temperature + z_gc_fraction + z_log_proteome_aa"
    original = fit_nb(data, formula)
    changed = fit_nb(data.assign(n_screened=data.n_screened * 2), formula)
    assert changed.params["Intercept"] == pytest.approx(original.params["Intercept"] - np.log(2), abs=1e-4)
    assert changed.params["z_trait_temperature"] == pytest.approx(original.params["z_trait_temperature"], abs=1e-4)


def test_phylogenetic_likelihood_builds_and_has_finite_logp(bundle):
    pytest.importorskip("pymc")
    tables, spec, _, tree = bundle
    cohort, _, _ = prepare_cohort(tables, spec)
    data, _ = analysis_subset(cohort)
    data = data.iloc[:8]
    design = pd.DataFrame({"Intercept": np.ones(len(data)), "x": np.linspace(-1, 1, len(data))})
    covariance = tree_covariance(tree, data.genome_id)
    model = build_phylogenetic_nb(data, design, covariance)
    assert np.isfinite(model.compile_logp()(model.initial_point()))


def test_observed_amrfinder_schema_is_normalized_without_guessing():
    import importlib.util
    module_spec = importlib.util.spec_from_file_location("amr_adapter", ROOT / "scripts/prepare_amr_subset.py")
    adapter = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(adapter)
    raw = pd.DataFrame({"Name": ["a", "b"], "Protein id": ["-", "protein2"], "Contig id": ["contig1", "contig2"], "Start": [1, 10], "Stop": [9, 90], "Strand": ["+", "-"], "Element symbol": ["gene1", "gene2"], "Type": ["AMR", "STRESS"], "Class": ["class1", "class2"], "Scope": ["core", "plus"]})
    result = adapter.normalize_chunk(raw, {"a", "b"}, core_only=True)
    assert result.genome_id.tolist() == ["a"]
    assert result.locus_id.tolist() == ["contig1:1:9:+"]
    assert result.gene_symbol.tolist() == ["gene1"]
    with pytest.raises(ValueError, match="missing columns"):
        adapter.normalize_chunk(raw.rename(columns={"Type": "Element type"}), {"a"})
