"""Count a released Bakta proteome's windows; never claim original AMP membership."""
import argparse
from datetime import datetime,timezone
import gzip
from importlib.metadata import version
import json
from pathlib import Path

import pandas as pd
from genomespot_lens.archaeome_archive import distinct_fragment_count
from genomespot_lens.core import sha256


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--annotation',type=Path,required=True)
    p.add_argument('--source-provenance',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    source=json.loads(a.source_provenance.read_text())
    if sha256(a.annotation)!=source['member_sha256']:raise ValueError('Annotation checksum mismatch')
    if not source['archive_md5_verified']:raise ValueError('Unverified source archive')
    doc=json.loads(a.annotation.read_text());proteins=[]
    for feature in doc['features']:
        if feature.get('aa'):
            if feature['type'] not in ['cds','sorf']:raise ValueError('Unexpected amino-acid feature type')
            proteins.append((feature['locus'],feature['aa']))
    if not proteins or len({p[0] for p in proteins})!=len(proteins):raise ValueError('Invalid protein registry')
    payload=''.join(f'>{name}\n{sequence}\n' for name,sequence in proteins).encode('ascii')
    (a.out/'released-proteins.fasta.gz').write_bytes(gzip.compress(payload,mtime=0))
    pd.DataFrame([{'protein_id':n,'length':len(s)} for n,s in proteins]).to_csv(a.out/'protein-registry.csv',index=False)
    sequences=[s for _,s in proteins]
    counts={str(n):distinct_fragment_count(sequences,8,n) for n in [49,50]}
    audit={'data_kind':'real','created_utc':datetime.now(timezone.utc).isoformat(),
        'source':source,'bakta_version':doc.get('version'),
        'sample_id':Path(source['selected_member']).name.removesuffix('.bakta.json'),'protein_records':len(sequences),
        'unique_protein_sequences':len(set(sequences)),'total_amino_acids':sum(map(len,sequences)),
        'distinct_canonical_fragments_min8_max':counts,'paper_reported_fragments':3919096,
        'single_sample_to_reported_ratio_max50':counts['50']/3919096,
        'original_screen_membership_known':False,'is_original_screen_export':False,
        'inference_allowed':'Conditional scope consistency only; no AMP scores or zero-hit genomes inferred',
        'analysis_code_sha256':sha256(Path(__file__)),
        'counter_code_sha256':sha256(Path(__file__).parents[1]/'src/genomespot_lens/archaeome_archive.py'),
        'packages':{n:version(n) for n in ['numpy','pydivsufsort']},
        'outputs':{f.name:sha256(f) for f in a.out.glob('*') if f.is_file() and f.name!='audit.json'}}
    (a.out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({k:audit[k] for k in ['sample_id','protein_records','distinct_canonical_fragments_min8_max','single_sample_to_reported_ratio_max50','original_screen_membership_known']},indent=2))


if __name__=='__main__':main()
