"""Independent Macrel small-ORF ecology comparison; not APEX replication."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from importlib.metadata import version

import numpy as np
import pandas as pd

from genomespot_lens.core import sha256
from genomespot_lens.ecology_comparison import fit_count_contrasts, within_species_permutation

SOURCE = 'https://doi.org/10.1186/s40168-025-02145-3'
SUPPLEMENT = 'https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs40168-025-02145-3/MediaObjects/40168_2025_2145_MOESM1_ESM.zip'


def prepare(folder):
    def read(n): return pd.read_excel(folder/f'Supplementary Data {n}.xlsx',header=1)
    registry=read(3);all_hits=read(4);hits=read(5)
    rename={'genome id in this study':'genome_id','source':'habitat',
        'completeness(%)':'completeness','contamination(%)':'contamination',
        'genome size (bp)':'genome_bp','GC content':'gc_fraction','contig N50':'contig_n50',
        'original MAG  id':'original_id','assembly status':'assembly_status'}
    registry=registry.rename(columns=rename)
    if registry.genome_id.duplicated().any() or hits.Access.duplicated().any():
        raise ValueError('Nonunique registry or positive identifiers')
    if not set(hits['genome id']).issubset(registry.genome_id):
        raise ValueError('Positive genomes absent from screened registry')
    mag_hits=all_hits[all_hits['genome id'].isin(registry.genome_id)]
    if set(hits.Access)!=set(mag_hits.Access): raise ValueError('S4/S5 positive mismatch')
    left=hits.drop(columns='habitat').sort_values('Access').reset_index(drop=True)
    right=mag_hits.sort_values('Access').reset_index(drop=True)
    pd.testing.assert_frame_equal(left,right)
    joined=hits.merge(registry[['genome_id','species','habitat']],
        left_on='genome id',right_on='genome_id',validate='many_to_one',suffixes=('_hit','_registry'))
    if not (joined.species_hit==joined.species_registry).all(): raise ValueError('Species mismatch')
    if not (joined.habitat_hit==joined.habitat_registry).all(): raise ValueError('Habitat mismatch')
    registry['amp_count']=registry.genome_id.map(hits.groupby('genome id').size()).fillna(0).astype(int)
    registry['multi_habitat_species']=registry.species.map(registry.groupby('species').habitat.nunique()).gt(1)
    if registry.isna().any().any(): raise ValueError('Incomplete registry')
    if (registry.genome_bp<=0).any(): raise ValueError('Invalid genome size')
    # Exact parent IDs and count integrity, not peptide design or rescoring.
    evidence=hits[['Access','GCF id','genome id','species','habitat']].rename(columns={'Access':'hit_id','genome id':'genome_id','GCF id':'published_family'})
    return registry,evidence


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    g,hits=prepare(a.raw)
    if len(g)!=1964 or len(hits)!=1327: raise ValueError('Published totals not reproduced')
    g['count_model_eligible']=g.multi_habitat_species & g.groupby('species').amp_count.transform('sum').gt(0)
    g.to_csv(a.out/'genomes.csv',index=False);hits.to_csv(a.out/'hit-registry.csv',index=False)
    summary=g.groupby('habitat').agg(genomes=('genome_id','size'),species=('species','nunique'),
        candidates=('amp_count','sum'),positive_genomes=('amp_count',lambda x:int(x.gt(0).sum())),
        mean_genome_bp=('genome_bp','mean'),mean_completeness=('completeness','mean'),
        mean_contamination=('contamination','mean'))
    summary['candidates_per_genome']=summary.candidates/summary.genomes
    summary.to_csv(a.out/'habitat-summary.csv')
    overlap=g.groupby(['species','habitat']).size().unstack(fill_value=0)
    overlap.to_csv(a.out/'species-habitat-overlap.csv')
    common=g[g.multi_habitat_species].copy()
    informative=common[common.count_model_eligible].copy()
    results=pd.concat([
        fit_count_contrasts(g,label='crude_all'),
        fit_count_contrasts(informative,label='crude_common_support'),
        fit_count_contrasts(informative,species_effect=True,label='species_common_support'),
        fit_count_contrasts(informative,species_effect=True,quality=True,label='species_quality_common_support')],ignore_index=True)
    results.to_csv(a.out/'effects.csv',index=False)
    permutation=within_species_permutation(common)
    (a.out/'permutation.json').write_text(json.dumps(permutation,indent=2)+'\n')
    provenance={'data_kind':'real','created_utc':datetime.now(timezone.utc).isoformat(),
        'exploratory_protocol_commit':'07a490c','article':SOURCE,'supplement':SUPPLEMENT,
        'license':'CC BY 4.0 (article and supplementary material; see article rights)',
        'endpoint':'Macrel small-ORF AMP candidates; distinct from encrypted APEX windows',
        'zero_definition':'Screened MAG registry members absent from complete positive table',
        'genomes':len(g),'positive_genomes':int(g.amp_count.gt(0).sum()),'zero_genomes':int(g.amp_count.eq(0).sum()),
        'hits':len(hits),'species':g.species.nunique(),'multi_habitat_species':common.species.nunique(),
        'informative_common_support_genomes':len(informative),
        'informative_common_support_species':informative.species.nunique(),
        'all_zero_common_species_genomes_conditioned_out':len(common)-len(informative),
        'all_zero_common_species_conditioned_out':common.species.nunique()-informative.species.nunique(),
        'common_support_genomes':len(common),'single_habitat_genomes':int((~g.multi_habitat_species).sum()),
        'screened_smorfs_per_genome_available':False,
        'covariates':'Standardized log genome bp, completeness, contamination, GC fraction; each centered/scaled within the fitted cohort',
        'uncertainty':'Species cluster sandwich, finite sample correction, t df=clusters-1; exploratory unadjusted 95% intervals',
        'raw_sha256':{f.name:sha256(f) for f in sorted(a.raw.glob('*.xlsx'))},
        'code_sha256':{'runner':sha256(Path(__file__)),
            'statistics':sha256(Path(__file__).parents[1]/'src/genomespot_lens/ecology_comparison.py')},
        'packages':{n:version(n) for n in ['numpy','pandas','statsmodels','scipy','openpyxl']}}
    provenance['supplement_zip_sha256']=sha256(a.raw.parent.parent/'lactobacillaceae-supplement.zip')
    provenance['outputs']={f.name:sha256(f) for f in a.out.glob('*') if f.is_file() and f.name!='provenance.json'}
    (a.out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(summary.to_string());print(results.to_string(index=False));print(json.dumps(permutation,indent=2))


if __name__=='__main__': main()
