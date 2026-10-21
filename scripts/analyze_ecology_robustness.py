"""Run frozen exploratory pair, quality, prevalence and species-influence checks."""
import argparse
from datetime import datetime,timezone
from importlib.metadata import version
import json
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from genomespot_lens.core import sha256
from genomespot_lens.ecology_comparison import fit_count_contrasts
from genomespot_lens.ecology_robustness import fit_direct_pair,PAIRS,resolved_species


def run(data,out):
    out.mkdir(parents=True,exist_ok=True)
    source_provenance=json.loads((data/'provenance.json').read_text())
    if source_provenance['data_kind']!='real':raise ValueError('Real registry required')
    if sha256(data/'genomes.csv')!=source_provenance['outputs']['genomes.csv']:
        raise ValueError('Changed baseline registry')
    g=pd.read_csv(data/'genomes.csv');results=[];members=[]
    for high in [False,True]:
        for presence in [False,True]:
            family=[]
            for numerator,denominator in PAIRS:
                r=fit_direct_pair(g,numerator,denominator,high_quality=high,presence=presence)
                for genome_id in r.pop('genome_ids'):
                    members.append({'endpoint':r['endpoint'],'quality':r['quality'],
                                    'contrast':r['contrast'],'genome_id':genome_id})
                r['species_ids']=';'.join(r['species_ids']);family.append(r)
            valid=[r for r in family if r['status']=='fit_pass']
            # Keep a three-test family even if a contrast cannot be estimated.
            adjusted=multipletests([r['p_value'] for r in valid]+[1.]*(3-len(valid)),method='holm')[1]
            for r,p in zip(valid,adjusted):r['p_holm_three_contrasts']=float(p)
            results+=family
    effects=pd.DataFrame(results);effects.to_csv(out/'direct-effects.csv',index=False)
    pd.DataFrame(members).to_csv(out/'cohort-members.csv',index=False)
    g['species_resolved']=resolved_species(g.species)
    g['resolved_joint_model_eligible']=g.count_model_eligible & g.species_resolved
    g.to_csv(out/'registry-with-taxonomy-status.csv',index=False)
    g[~g.species_resolved].to_csv(out/'unresolved-species-exclusions.csv',index=False)
    matched=g[g.resolved_joint_model_eligible].copy()
    corrected=pd.concat([
        fit_count_contrasts(g[g.species_resolved],label='crude_resolved_registry'),
        fit_count_contrasts(matched,label='crude_resolved_common'),
        fit_count_contrasts(matched,species_effect=True,label='species_resolved_common'),
        fit_count_contrasts(matched,species_effect=True,quality=True,label='species_quality_resolved_common')],ignore_index=True)
    corrected.to_csv(out/'corrected-joint-effects.csv',index=False)
    influence=[]
    for species in sorted(matched.species.unique()):
        remaining=matched[matched.species!=species].copy()
        fit=fit_count_contrasts(remaining,species_effect=True,quality=True,label='leave_one_species_out')
        fit['excluded_species']=species
        fit['excluded_genomes']=len(matched)-len(remaining)
        influence.append(fit)
    influence=pd.concat(influence,ignore_index=True)
    influence.to_csv(out/'species-influence.csv',index=False)
    summary=influence.groupby('contrast').agg(fits=('rate_ratio','size'),
        minimum_rate_ratio=('rate_ratio','min'),maximum_rate_ratio=('rate_ratio','max'),
        minimum_ci_low=('ci_low','min'),maximum_ci_high=('ci_high','max'))
    summary['intervals_excluding_one']=influence.assign(excludes=(influence.ci_low>1)|(influence.ci_high<1)).groupby('contrast').excludes.sum()
    summary.to_csv(out/'influence-summary.csv')
    provenance={'data_kind':'real','created_utc':datetime.now(timezone.utc).isoformat(),
        'exploratory_protocol_commit':'d5c1eac','taxonomy_correction_protocol_commit':'78d5882','baseline_genomes_sha256':sha256(data/'genomes.csv'),
        'baseline_provenance_sha256':sha256(data/'provenance.json'),
        'software':{n:version(n) for n in ['numpy','pandas','scipy','statsmodels']},
        'analysis_code_sha256':sha256(Path(__file__)),
        'new_statistics_sha256':sha256(Path(__file__).parents[1]/'src/genomespot_lens/ecology_robustness.py'),
        'original_statistics_sha256':sha256(Path(__file__).parents[1]/'src/genomespot_lens/ecology_comparison.py'),
        'counts':{'registry_genomes':len(g),'unresolved_species_genomes':int((~g.species_resolved).sum()),'resolved_species':g[g.species_resolved].species.nunique(),'direct_specifications':len(results),
                  'informative_joint_species':matched.species.nunique(),'leave_species_out_fits':len(influence)//3},
        'quality_definition':'completeness >=90%; contamination <5%; overlap recomputed after filtering',
        'screened_smorfs_denominator_available':False,
        'inference':'Exploratory species-cluster t intervals; direct-pair within-species estimand; no cross-species phylogenetic correction',
        'outputs':{f.name:sha256(f) for f in out.glob('*.csv')}}
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(effects.drop(columns='species_ids').to_string(index=False));print(summary.to_string())
    return effects,influence


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.data,a.out)
