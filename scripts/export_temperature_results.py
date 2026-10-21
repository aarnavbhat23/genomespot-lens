"""Export only diagnostic-passing final fits, with a separate sampling-attempt ledger."""
import argparse
import json
from pathlib import Path
import shutil

import pandas as pd
from genomespot_lens.core import sha256


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline-run',type=Path,required=True);p.add_argument('--other-run',type=Path,required=True)
    p.add_argument('--earlier-run',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    results=[];sources={}
    for name in ['baseline','adjusted','parents100','broad','predicted']:
        source=a.baseline_run if name=='baseline' else a.other_run
        result=json.loads((source/f'{name}_summary.json').read_text())
        if result['status']!='diagnostics_pass' or result['divergences'] or result['max_rhat']>=1.01 or result['minimum_bulk_ess']<400:
            raise ValueError(f'Final model {name} has not passed the protocol diagnostics')
        results.append(result)
        for suffix in ['summary.json','diagnostics.csv','predictive.csv']:
            shutil.copy2(source/f'{name}_{suffix}',a.out/f'{name}_{suffix}')
        sources[name]={'posterior_draws_sha256':sha256(source/f'{name}.nc'),
                      'run_provenance':json.loads((source/'run_provenance.json').read_text())}
    (a.out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    pd.DataFrame(results).drop(columns=['aa_pca_variance_fraction','temperature_range',
        'posterior_predictive_zero_count_q025_q975','posterior_predictive_total_q025_q975']).to_csv(a.out/'effects.csv',index=False)
    earlier=json.loads((a.earlier_run/'results.json').read_text())
    (a.out/'sampling-attempts.json').write_text(json.dumps({'earlier_attempts':earlier,
        'final_source_records':sources,'selection_rule':'Same models; final stricter sampling configuration selected by protocol diagnostics, never significance. Earlier intervals with failed diagnostics are not reported as findings.'},indent=2)+'\n')
    (a.out/'provenance.json').write_text(json.dumps({'data_kind':'real','all_models_passed_diagnostics':True,
        'posterior_draws_public':False,'posterior_draws_location':'Local runs/temperature-repair-*/*.nc',
        'export_code_sha256':sha256(Path(__file__)),
        'shared_model_code_sha256':sha256(Path(__file__).resolve().parents[1]/'src/genomespot_lens/amp_study.py'),
        'source_runs':sources,'outputs':{p.name:sha256(p) for p in a.out.iterdir()}},indent=2)+'\n')
    print('Exported five diagnostic-passing real-data models')


if __name__=='__main__':main()
