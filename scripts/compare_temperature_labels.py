"""Compare labels and disclose overlap; never call this independent validation."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from genomespot_lens.core import sha256


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True);p.add_argument('--upstream',type=Path,required=True)
    a=p.parse_args();g=pd.read_csv(a.data/'genomes.csv');phenotypes=pd.read_csv(a.data/'phenotypes.csv')
    d=g[g.primary_eligible & g.predicted_temperature_unflagged].merge(phenotypes[['genome_id','matched_assembly_ids']],on='genome_id')
    base=a.upstream/'data/holdouts'
    files=[base/'train_set_temperature.txt',base/'test_set_temperature.txt']
    for name,path in zip(['listed_training_accession_overlap','listed_original_holdout_accession_overlap'],files):
        listed=set(path.read_text().split())
        # The source lists omit versions. This flags potential overlap only;
        # stripping versions is never allowed in the source/phenotype mapping.
        d[name]=d.matched_assembly_ids.map(lambda s:bool({x.split('.')[0] for x in s.split(';')}&listed))
    d[['genome_id','matched_assembly_ids','listed_training_accession_overlap','listed_original_holdout_accession_overlap']].to_csv(a.data/'training-overlap.csv',index=False)
    result={'genomes':len(d),'mae_celsius':float(np.mean(abs(d.predicted_temperature-d.observed_temperature))),
        'pearson_correlation':float(d.predicted_temperature.corr(d.observed_temperature)),
        'listed_training_or_original_holdout_overlap':int((d.listed_training_accession_overlap|d.listed_original_holdout_accession_overlap).sum()),
        'interpretation':'Label agreement only, NOT independent predictive validation. Final GenomeSPOT models use both original training and holdout sets; accession lists omit versions, so overlap audit is conservative and never used for phenotype joining.',
        'holdout_file_sha256':{path.name:sha256(path) for path in files}}
    (a.data/'temperature-agreement.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
