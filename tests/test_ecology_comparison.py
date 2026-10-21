"""Check real registry integrity and independent properties of species controls."""
import importlib.util
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
import pytest

pytest.importorskip('statsmodels')
pytest.importorskip('openpyxl')
from genomespot_lens.core import sha256
from genomespot_lens.ecology_comparison import within_species_permutation

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'examples/lactobacillaceae-ecology'


def load_script(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_full_source_registry_reproduces_zero_and_positive_genomes(tmp_path):
    zip_path=DATA/'sources/lactobacillaceae-supplement.zip'
    provenance=json.loads((DATA/'provenance.json').read_text())
    assert sha256(zip_path)==provenance['supplement_zip_sha256']
    with zipfile.ZipFile(zip_path) as archive: archive.extractall(tmp_path)
    folder=tmp_path/'Additional File 1'
    for name,digest in provenance['raw_sha256'].items(): assert sha256(folder/name)==digest
    g,hits=load_script('analyze_lactobacillaceae_ecology').prepare(folder)
    public=pd.read_csv(DATA/'genomes.csv')
    assert len(g)==1964 and len(hits)==1327
    assert g.amp_count.eq(0).sum()==1084
    assert g.amp_count.gt(0).sum()==880
    assert g.amp_count.sum()==1327
    pd.testing.assert_series_equal(g.amp_count,public.amp_count)
    pd.testing.assert_frame_equal(hits,pd.read_csv(DATA/'hit-registry.csv'))
    assert public.multi_habitat_species.sum()==1323
    assert public[public.multi_habitat_species].species.nunique()==40
    assert public.count_model_eligible.sum()==1312
    assert public[public.count_model_eligible].species.nunique()==36


def test_exported_comparison_inputs_and_matched_cohorts_are_frozen():
    p=json.loads((DATA/'provenance.json').read_text())
    for name,digest in p['outputs'].items(): assert sha256(DATA/name)==digest
    assert sha256(ROOT/'scripts/analyze_lactobacillaceae_ecology.py')==p['code_sha256']['runner']
    assert sha256(ROOT/'src/genomespot_lens/ecology_comparison.py')==p['code_sha256']['statistics']
    fits=pd.read_csv(DATA/'effects.csv')
    matched=fits[fits.model!='crude_all']
    assert set(matched.genomes)=={1312}
    assert set(matched.species_clusters)=={36}
    assert fits.converged.all() and np.isfinite(fits.select_dtypes('number')).all().all()
    assert (fits.ci_low<fits.rate_ratio).all() and (fits.rate_ratio<fits.ci_high).all()


def test_conditional_score_matches_hand_calculation_and_exact_null_probability():
    # Each species has one genome/habitat and counts 0,0,6. The score is 24;
    # of 3^2 equiprobable placements of the two nonzero counts, 3 tie the max.
    toy=pd.DataFrame({'species':['a']*3+['b']*3,
        'habitat':['animal','human','fermented food']*2,'amp_count':[0,0,6]*2})
    result=within_species_permutation(toy,permutations=2999,seed=23)
    assert result['statistic']==pytest.approx(24)
    assert result['p_value']==pytest.approx(1/3,abs=.04)
    zero=toy.assign(amp_count=2)
    assert within_species_permutation(zero,permutations=19)['p_value']==1


def test_habitat_confounding_without_within_species_overlap_is_not_identifiable():
    toy=pd.DataFrame({'species':['a']*3+['b']*3,
        'habitat':['animal']*3+['human']*3,'amp_count':[1]*3+[10]*3})
    with pytest.raises(ValueError,match='No full habitat contrast'):
        within_species_permutation(toy,permutations=19)


def test_real_permutation_reproduces_counts_with_score_roundoff_tolerance():
    g=pd.read_csv(DATA/'genomes.csv')
    observed=within_species_permutation(g[g.multi_habitat_species],permutations=9999,seed=23)
    expected=json.loads((DATA/'permutation.json').read_text())
    assert observed['statistic']==pytest.approx(expected['statistic'],rel=1e-12,abs=1e-12)
    assert {k:v for k,v in observed.items() if k!='statistic'}=={k:v for k,v in expected.items() if k!='statistic'}
