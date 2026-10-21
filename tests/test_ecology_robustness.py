"""Guard taxonomy resolution, direct support, endpoint identity and real artifacts."""
from pathlib import Path
import gzip,json
import numpy as np
import pandas as pd
import pytest
pytest.importorskip('statsmodels')
pytest.importorskip('pydivsufsort')
from Bio import SeqIO
from genomespot_lens.core import sha256
from genomespot_lens.ecology_robustness import resolved_species,pair_cohort,fit_direct_pair
from genomespot_lens.archaeome_archive import distinct_fragment_count
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'examples/ecology-robustness'


def test_bare_species_rank_cannot_establish_species_membership():
    labels=pd.Series(['s__','s__  ',None,'s__unknown','s__Lactobacillus sp945980025','s__Pediococcus pentosaceus'])
    assert resolved_species(labels).tolist()==[False,False,False,False,True,True]


def test_overlap_is_recomputed_after_quality_and_individual_zeros_remain():
    g=pd.DataFrame({'genome_id':['a','b','c','d','e','f'],
      'species':['s__A one','s__A one','s__B two','s__B two','s__','s__'],
      'habitat':['human','animal']*3,'amp_count':[0,1,0,2,9,10],
      'completeness':[95,80,95,95,95,95],'contamination':[0]*6})
    full=pair_cohort(g,'human','animal');high=pair_cohort(g,'human','animal',high_quality=True)
    assert set(full.genome_id)=={'a','b','c','d'}
    assert set(high.genome_id)=={'c','d'}
    assert high.amp_count.eq(0).sum()==1


def test_presence_endpoint_does_not_reuse_constant_binary_species():
    g=pd.DataFrame({'genome_id':['a','b','c','d'],'species':['s__A one']*2+['s__B two']*2,
      'habitat':['human','animal']*2,'amp_count':[1,2,0,2]})
    assert set(pair_cohort(g,'human','animal').genome_id)=={'a','b','c','d'}
    assert set(pair_cohort(g,'human','animal',presence=True).genome_id)=={'c','d'}


def test_new_artifacts_pin_sources_and_analysis_code():
    p=json.loads((DATA/'provenance.json').read_text())
    for name,digest in p['outputs'].items():assert sha256(DATA/name)==digest
    assert sha256(ROOT/'scripts/analyze_ecology_robustness.py')==p['analysis_code_sha256']
    assert sha256(ROOT/'src/genomespot_lens/ecology_robustness.py')==p['new_statistics_sha256']
    g=pd.read_csv(DATA/'registry-with-taxonomy-status.csv')
    missing=g[~g.species_resolved]
    assert len(g)==1964 and len(missing)==20 and missing.genus.nunique()==6
    assert not missing.resolved_joint_model_eligible.any()
    matched=g[g.resolved_joint_model_eligible]
    assert len(matched)==1292 and matched.species.nunique()==35
    p=json.loads((DATA/'pair-influence/provenance.json').read_text())
    assert sha256(ROOT/'scripts/audit_ecology_pair_influence.py')==p['analysis_code_sha256']
    assert sha256(ROOT/'src/genomespot_lens/ecology_robustness.py')==p['statistics_sha256']
    for name,digest in p['outputs'].items():assert sha256(DATA/'pair-influence'/name)==digest


def test_all_direct_cohorts_have_real_resolved_species_in_both_habitats():
    g=pd.read_csv(ROOT/'examples/lactobacillaceae-ecology/genomes.csv')
    effects=pd.read_csv(DATA/'direct-effects.csv');members=pd.read_csv(DATA/'cohort-members.csv')
    assert len(effects)==12
    for row in effects.itertuples():
        numerator,denominator=row.contrast.split(' / ')
        expected=pair_cohort(g,numerator,denominator,high_quality=row.quality=='high',presence=row.endpoint=='presence')
        selected=members[(members.contrast==row.contrast)&(members.quality==row.quality)&(members.endpoint==row.endpoint)]
        assert set(selected.genome_id)==set(expected.genome_id)
        assert len(expected)==row.genomes and expected.species.nunique()==row.species_clusters
        assert resolved_species(expected.species).all()
        assert expected.groupby('species').habitat.nunique().eq(2).all()
        if row.status=='fit_pass':assert row.species_clusters>=10
        if row.status=='insufficient_support':assert row.species_clusters<10


def test_real_direct_effect_is_order_invariant_and_reproducible():
    g=pd.read_csv(ROOT/'examples/lactobacillaceae-ecology/genomes.csv')
    actual=fit_direct_pair(g.sample(frac=1,random_state=17),'human','animal',high_quality=True)
    expected=pd.read_csv(DATA/'direct-effects.csv')
    expected=expected[(expected.contrast=='human / animal')&(expected.quality=='high')&(expected.endpoint=='count')].iloc[0]
    assert actual['status']=='fit_pass' and actual['genomes']==379 and actual['species_clusters']==16
    for key in ['effect','ci_low','ci_high','p_value']:assert actual[key]==pytest.approx(expected[key],rel=1e-7,abs=1e-10)


def test_insufficient_support_is_not_a_zero_or_significant_effect():
    d=pd.read_csv(DATA/'direct-effects.csv')
    insufficient=d[d.status=='insufficient_support']
    assert len(insufficient)==2 and insufficient.effect.isna().all()
    withheld=d[d.status.str.startswith('withhold')]
    assert len(withheld)==1 and withheld.effect.isna().all()


def test_released_bakta_count_is_not_mislabeled_an_original_screen():
    data=ROOT/'examples/atb-screen-audit';audit=json.loads((data/'audit.json').read_text())
    assert audit['is_original_screen_export'] is False
    assert audit['original_screen_membership_known'] is False
    assert audit['analysis_code_sha256']==sha256(ROOT/'scripts/audit_atb_screen_scope.py')
    for name,digest in audit['outputs'].items():assert sha256(data/name)==digest
    with gzip.open(data/'released-proteins.fasta.gz','rt') as f:sequences=[str(r.seq) for r in SeqIO.parse(f,'fasta')]
    assert len(sequences)==audit['protein_records']==2773
    assert distinct_fragment_count(sequences)==audit['distinct_canonical_fragments_min8_max']['50']==33330305
    assert distinct_fragment_count(sequences,8,49)==audit['distinct_canonical_fragments_min8_max']['49']==32612241


def test_export_search_inventory_preserves_limits_and_frozen_counts():
    data=ROOT/'examples/atb-screen-audit'
    manifest=json.loads((data/'search-manifest.json').read_text())
    assert manifest['original_export_recovered'] is False
    assert manifest['limits']
    for name,digest in manifest['inventory_sha256'].items():assert sha256(data/name)==digest
    assert sha256(data/'osf-file-index.tsv')==manifest['osf_index_sha256']
    assert len(pd.read_csv(data/'osf-file-index.tsv',sep='\t'))==manifest['counts']['osf_indexed_files']==3009
    assert len(pd.read_csv(data/'gitlab-file-inventory.csv'))==manifest['counts']['gitlab_blobs_default_and_two_alternate_apex_branches']==933
    assert len(pd.read_csv(data/'huggingface-repositories.csv'))==10
