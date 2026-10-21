"""Post-lead influence audit of every direct-pair count specification."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import pandas as pd
from genomespot_lens.core import sha256
from genomespot_lens.ecology_robustness import fit_direct_pair,pair_cohort,PAIRS


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    g=pd.read_csv(a.data/'genomes.csv');rows=[];support=[]
    for high in [False,True]:
        for numerator,denominator in PAIRS:
            cohort=pair_cohort(g,numerator,denominator,high_quality=high)
            table=cohort.groupby(['species','habitat']).agg(genomes=('genome_id','size'),
                candidates=('amp_count','sum'),positive_genomes=('amp_count',lambda x:int(x.gt(0).sum())))
            table['mean_candidates']=table.candidates/table.genomes
            table=table.reset_index();table['quality']='high' if high else 'published'
            table['contrast']=f'{numerator} / {denominator}';support.append(table)
            for excluded in sorted(cohort.species.unique()):
                r=fit_direct_pair(g[g.species!=excluded],numerator,denominator,high_quality=high)
                r.pop('genome_ids');r.pop('species_ids');r['excluded_species']=excluded
                rows.append(r)
    pd.concat(support,ignore_index=True).to_csv(a.out/'species-habitat-support.csv',index=False)
    results=pd.DataFrame(rows);results.to_csv(a.out/'direct-pair-influence.csv',index=False)
    valid=results[results.status=='fit_pass'].copy()
    summaries=[]
    for (quality,contrast),d in results.groupby(['quality','contrast']):
        ok=d[d.status=='fit_pass']
        summaries.append({'quality':quality,'contrast':contrast,'omissions':len(d),
            'estimable_omissions':len(ok),'insufficient_support':int(d.status.eq('insufficient_support').sum()),
            'minimum_effect':ok.effect.min() if len(ok) else None,
            'maximum_effect':ok.effect.max() if len(ok) else None,
            'intervals_below_one':int(ok.ci_high.lt(1).sum()),
            'intervals_above_one':int(ok.ci_low.gt(1).sum())})
    summary=pd.DataFrame(summaries);summary.to_csv(a.out/'pair-influence-summary.csv',index=False)
    provenance={'data_kind':'real','created_utc':datetime.now(timezone.utc).isoformat(),
        'baseline_genomes_sha256':sha256(a.data/'genomes.csv'),'protocol_commit':'405ac97',
        'taxonomy_correction_protocol_commit':'78d5882','post_lead_analysis':True,'analysis_code_sha256':sha256(Path(__file__)),
        'statistics_sha256':sha256(Path(__file__).parents[1]/'src/genomespot_lens/ecology_robustness.py'),
        'outputs':{f.name:sha256(f) for f in a.out.glob('*.csv')}}
    (a.out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(summary.to_string(index=False))


if __name__=='__main__':main()
