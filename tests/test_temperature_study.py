"""Verify strain joins and exact subset accounting, including real-data artifacts."""
import importlib.util
from pathlib import Path

import pytest
pytest.importorskip('ahocorasick')
pytest.importorskip('pydivsufsort')

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('temperature_preparation',ROOT/'scripts/prepare_temperature_study.py')
preparation=importlib.util.module_from_spec(spec);spec.loader.exec_module(preparation)


def test_culture_aliases_do_not_turn_species_or_assembly_into_strain_identity():
    rows=[{'strain':'8Ac'}]
    ncbi={'assembly_info':{'biosample':{'strain':'6Ac'}}}
    assert not preparation.strain_supported(rows,'Methanothrix harundinacea (strain 6Ac).',ncbi)
    assert preparation.strain_supported([{'strain':'∆H'}], 'Methanothermobacter (strain Delta H).',{})
    assert preparation.strain_supported([{'strain':'Z'}], 'Methanocorpusculum (strain DSM 4855 / Z).',{})
    assert not preparation.strain_supported([{'strain':'Hvv3'}], 'Thermofilum (strain Hrk 5).',{})
    assert preparation.strain_supported([{'strain':'Hvv3'}], 'Thermofilum (strain Hrk 5).',{},['Hvv3'])
    assert not preparation.strain_supported([], 'any', {})


def test_jsonp_is_data_and_complete_page_is_required(tmp_path):
    import json
    row={'id':1,'cell':['-']*21}
    row['cell'][2]='<a href="/taxon/123">123</a>'
    path=tmp_path/'tempura.json';path.write_text('/**/callback('+json.dumps({'total':1,'rows':[row]})+')')
    parsed=preparation.parse_tempura(path)
    assert parsed[0]['taxid']=='123'
    path.write_text('/**/callback('+json.dumps({'total':2,'rows':[row]})+')')
    with pytest.raises(ValueError,match='Incomplete'):preparation.parse_tempura(path)


def test_actual_subset_candidate_unions_do_not_transfer_unmatched_parent_hits():
    import csv,json,hashlib
    from collections import defaultdict
    data=ROOT/'examples/temperature-study'
    provenance=json.loads((data/'provenance.json').read_text())
    for name in ['genomes.csv','proteins.csv','phenotypes.csv','exclusions.csv','phylogeny.nwk']:
        assert hashlib.sha256((data/name).read_bytes()).hexdigest()==provenance['outputs'][name]
    genomes=list(csv.DictReader((data/'genomes.csv').open()))
    proteins=list(csv.DictReader((data/'proteins.csv').open()))
    assert len(proteins)==15786
    assert sum(x['primary_eligible']=='True' for x in genomes)==29
    assert len({(r['taxon_id'],r['sequence_sha256']) for r in proteins})==len(proteins)
    source=ROOT/'examples/amp-study/archaeome-reconstruction/positive_parent_mappings.csv'
    positive=defaultdict(set)
    for row in csv.DictReader(source.open()):positive[(row['taxid'],row['accession'])].add(row['catalog_row'])
    unions=defaultdict(set);exposures=defaultdict(int)
    for p in proteins:
        hits=positive[(p['taxon_id'],p['parent_id'])]
        assert len(hits)==int(p['amp_count'])
        unions[p['genome_id']].update(hits);exposures[p['genome_id']]+=int(p['n_screened'])
    for g in genomes:
        assert len(unions[g['genome_id']])==int(g['amp_count'])
        assert int(g['n_screened'])<=exposures[g['genome_id']]
        if g['primary_eligible']=='True':
            assert g['temperature_strain_supported']=='True'
            assert g['genome_qc_pass']=='True'
            assert int(g['matched_source_proteins'])>=20
            assert g['gtdb_tip']


def test_pruned_marker_tree_preserves_original_patristic_distances():
    import csv,gzip
    from Bio import Phylo
    data=ROOT/'examples/temperature-study'
    with gzip.open(data/'sources/gtdb-ar53-r226.tree.gz','rt') as f:original=Phylo.read(f,'newick')
    pruned=Phylo.read(data/'phylogeny.nwk','newick')
    mapped={r['genome_id']:r['gtdb_tip'] for r in csv.DictReader((data/'genomes.csv').open()) if r['gtdb_tip']}
    assert {t.name for t in pruned.get_terminals()}==set(mapped)
    ids=sorted(mapped)[::11]
    for a,b in zip(ids,ids[1:]):
        assert pruned.distance(a,b)==pytest.approx(original.distance(mapped[a],mapped[b]),abs=1e-8)


def test_posterior_hdf_backend_round_trip(tmp_path):
    az=pytest.importorskip('arviz')
    import numpy as np
    sample=az.from_dict(posterior={'beta':np.arange(20).reshape(2,10)})
    path=tmp_path/'posterior.nc';sample.to_netcdf(path)
    assert np.array_equal(az.from_netcdf(path).posterior.beta.values,sample.posterior.beta.values)


def test_exported_real_posterior_claims_pass_the_original_diagnostic_gates():
    import json,hashlib
    folder=ROOT/'examples/temperature-study/posterior-results'
    results=json.loads((folder/'results.json').read_text())
    assert {r['model'] for r in results}=={'baseline','adjusted','parents100','broad','predicted'}
    for r in results:
        assert r['status']=='diagnostics_pass'
        assert r['max_rhat']<1.01 and r['minimum_bulk_ess']>=400 and r['divergences']==0
        assert 0<r['rr_q025']<r['rr_per10C_median']<r['rr_q975']
        assert r['draws_per_chain']==4000 and r['tune_per_chain']==3000
    provenance=json.loads((folder/'provenance.json').read_text())
    current=hashlib.sha256((ROOT/'scripts/analyze_temperature_study.py').read_bytes()).hexdigest()
    data_hash=hashlib.sha256((ROOT/'examples/temperature-study/provenance.json').read_bytes()).hexdigest()
    for source in provenance['source_runs'].values():
        assert source['run_provenance']['analysis_code_sha256']==current
        assert source['run_provenance']['data_provenance_sha256']==data_hash
    by={r['model']:r for r in results}
    assert by['baseline']['genomes']==by['adjusted']['genomes']==29


def test_tree_release_sensitivity_preserves_versioned_tips_and_distances(monkeypatch):
    import gzip,json
    import pandas as pd
    from Bio import Phylo
    monkeypatch.syspath_prepend(str(ROOT/'scripts'))
    import analyze_temperature_sensitivity as sensitivity
    data=ROOT/'examples/temperature-study'
    source=data/'sources/ar53-r232.tree.gz'
    with gzip.open(source,'rt') as f:original=Phylo.read(f,'newick')
    pruned=Phylo.read(data/'sensitivity-results/phylogeny-r232.nwk','newick')
    actual={t.name for t in pruned.get_terminals()}
    p=json.loads((data/'sensitivity-results/provenance.json').read_text())
    assert actual==set(p['primary_genome_ids'])==set(p['r232_common_genome_ids'])
    mapping={t.name.removeprefix('RS_').removeprefix('GB_'):t.name for t in original.get_terminals()}
    names=sorted(actual)[::4]
    for a,b in zip(names,names[1:]):
        assert pruned.distance(a,b)==pytest.approx(original.distance(mapping[a],mapping[b]),abs=1e-8)


def test_new_sensitivity_claims_pass_diagnostics_without_replacing_originals():
    import json,hashlib
    data=ROOT/'examples/temperature-study'
    p=json.loads((data/'sensitivity-results/provenance.json').read_text())
    digest=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
    assert p['analysis_code_sha256']==digest(ROOT/'scripts/analyze_temperature_sensitivity.py')
    assert p['original_design_code_sha256']==digest(ROOT/'scripts/analyze_temperature_study.py')
    assert p['r232_source_sha256']==digest(data/'sources/ar53-r232.tree.gz')
    assert p['r232_pruned_sha256']==digest(data/'sensitivity-results/phylogeny-r232.nwk')
    for name,value in p['data_sha256'].items():assert digest(data/name)==value
    results=json.loads((data/'sensitivity-results/results.json').read_text())
    assert {r['model'] for r in results}=={'prior075','prior300','tree232'}
    for r in results:
        assert r['status']=='diagnostics_pass' and r['divergences']==0
        assert r['max_rhat']<1.01 and r['minimum_bulk_ess']>=400
        assert r['genomes']==29 and set(r['genome_ids'])==set(p['primary_genome_ids'])


def test_published_root_pruning_retains_shared_stems_not_just_pair_distances(monkeypatch):
    import gzip,json
    from Bio import Phylo
    monkeypatch.syspath_prepend(str(ROOT/'scripts'))
    import analyze_temperature_rooting
    data=ROOT/'examples/temperature-study'
    p=json.loads((data/'rooting-results/provenance.json').read_text())
    for label,source in [('r226_published_root','gtdb-ar53-r226.tree.gz'),('r232_published_root','ar53-r232.tree.gz')]:
        with gzip.open(data/'sources'/source,'rt') as f:original=Phylo.read(f,'newick')
        pruned=Phylo.read(data/'rooting-results'/f'{label}.nwk','newick')
        assert len(pruned.root.clades)==1
        mapping={t.name.removeprefix('RS_').removeprefix('GB_'):t.name for t in original.get_terminals()}
        names=p['genome_ids'];raw_names=[mapping[n] for n in names]
        stem=original.distance(original.root,original.common_ancestor(raw_names))
        assert stem>0
        assert pruned.distance(pruned.root,pruned.common_ancestor(names))==pytest.approx(stem,abs=1e-8)
        for name in names[::5]:
            assert pruned.distance(pruned.root,name)==pytest.approx(original.distance(original.root,mapping[name]),abs=1e-8)


def test_published_root_fits_are_diagnostic_gated_and_source_pinned():
    import json,hashlib
    data=ROOT/'examples/temperature-study'
    folder=data/'rooting-results'
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    p=json.loads((folder/'provenance.json').read_text())
    assert p['analysis_code_sha256']==digest(ROOT/'scripts/analyze_temperature_rooting.py')
    assert p['shared_model_code_sha256']==digest(ROOT/'scripts/analyze_temperature_sensitivity.py')
    assert p['original_design_code_sha256']==digest(ROOT/'scripts/analyze_temperature_study.py')
    assert p['genomes_sha256']==digest(data/'genomes.csv')
    for label,raw in [('r226_published_root','gtdb-ar53-r226.tree.gz'),('r232_published_root','ar53-r232.tree.gz')]:
        assert p['raw_tree_sha256'][label]==digest(data/'sources'/raw)
    results=json.loads((folder/'results.json').read_text())
    assert {r['model'] for r in results}=={'r226_published_root','r232_published_root'}
    for r in results:
        assert r['status']=='diagnostics_pass' and r['genomes']==29
        assert r['divergences']==0 and r['max_rhat']<1.01 and r['minimum_bulk_ess']>=400
        assert r['tree_sha256']==digest(folder/f"{r['model']}.nwk")
