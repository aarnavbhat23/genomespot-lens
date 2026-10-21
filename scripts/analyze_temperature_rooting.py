"""Retain published roots after detecting shared-stem removal in earlier pruning."""
import argparse
from datetime import datetime, timezone
import gzip
import json
from pathlib import Path

import arviz as az
from Bio import Phylo
import numpy as np
import pandas as pd
import pymc as pm

from analyze_temperature_study import design
from analyze_temperature_sensitivity import sensitivity_model
from genomespot_lens.amp_study import tree_covariance
from genomespot_lens.core import sha256


def prune_published_root(source, destination, genome_ids):
    with gzip.open(source,'rt') as f: tree=Phylo.read(f,'newick')
    selected=set(genome_ids)
    aliases={t.name:t.name.removeprefix('RS_').removeprefix('GB_') for t in tree.get_terminals()}
    mapping={tip:asm for tip,asm in aliases.items() if asm in selected}
    if len(mapping)!=len(selected) or len(set(mapping.values()))!=len(selected):
        raise ValueError('Missing or ambiguous exact versioned assembly tips')
    depths={asm:tree.distance(tree.root,tip) for tip,asm in mapping.items()}
    shared_stem=tree.distance(tree.root,tree.common_ancestor(list(mapping)))
    def retain(node,is_root=False):
        if node.is_terminal():
            if node.name not in mapping:return None
            node.name=mapping[node.name];return node
        node.clades=[child for c in node.clades if (child:=retain(c)) is not None]
        if not node.clades:return None
        if len(node.clades)==1 and not is_root:
            child=node.clades[0]
            child.branch_length=(child.branch_length or 0)+(node.branch_length or 0)
            return child
        return node
    tree.root=retain(tree.root,is_root=True);tree.root.branch_length=0;tree.rooted=True
    Phylo.write(tree,destination,'newick',format_branch_length='%1.10f')
    recovered=Phylo.read(destination,'newick')
    if not all(np.isclose(recovered.distance(recovered.root,k),v,rtol=0,atol=1e-8) for k,v in depths.items()):
        raise ValueError('Root-to-tip depth was changed')
    if not np.isclose(recovered.distance(recovered.root,recovered.common_ancestor(sorted(selected))),shared_stem,rtol=0,atol=1e-8):
        raise ValueError('Shared root stem was changed')
    return float(shared_stem)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    data=pd.read_csv(a.data/'genomes.csv');data=data[data.primary_eligible].reset_index(drop=True)
    if len(data)<20:raise ValueError('Below feasibility threshold')
    x,design_info=design(data,'observed_temperature',True)
    sources={'r226_published_root':a.data/'sources/gtdb-ar53-r226.tree.gz',
             'r232_published_root':a.data/'sources/ar53-r232.tree.gz'}
    provenance={'created_utc':datetime.now(timezone.utc).isoformat(),'data_kind':'real',
        'protocol_commit':'b089cb9','analysis_code_sha256':sha256(Path(__file__)),
        'shared_model_code_sha256':sha256(Path(__file__).with_name('analyze_temperature_sensitivity.py')),
        'original_design_code_sha256':sha256(Path(__file__).with_name('analyze_temperature_study.py')),
        'genomes_sha256':sha256(a.data/'genomes.csv'),'genome_ids':data.genome_id.tolist(),
        'raw_tree_sha256':{k:sha256(v) for k,v in sources.items()},
        'rooting':'Full published root retained, including shared unary stem; no independent rerooting',
        'pymc_version':pm.__version__,'arviz_version':az.__version__}
    results=[]
    for label,source in sources.items():
        tree=a.out/f'{label}.nwk'
        stem=prune_published_root(source,tree,data.genome_id)
        model=sensitivity_model(data,x,tree_covariance(tree,data.genome_id.tolist()),1.5)
        print('Fitting',label,'shared stem',stem,flush=True)
        with model:
            fit=pm.sample(draws=4000,tune=3000,chains=4,cores=2,target_accept=.995,
                random_seed=23,progressbar=False,idata_kwargs={'log_likelihood':True})
        fit.to_netcdf(a.out/f'{label}.nc')
        az.summary(fit,var_names=['beta','phylo_sd','nb_shape']).to_csv(a.out/f'{label}_diagnostics.csv')
        rhat=float(az.rhat(fit).to_array().max());ess=float(az.ess(fit,method='bulk').to_array().min())
        div=int(fit.sample_stats.diverging.sum())
        beta=fit.posterior.beta.sel(coefficient='temperature_per10C').values.ravel()
        ci=np.quantile(np.exp(beta),[.025,.5,.975])
        result={'model':label,'genomes':len(data),'shared_stem_substitutions_per_site':stem,
            'rr_q025':float(ci[0]),'rr_per10C_median':float(ci[1]),'rr_q975':float(ci[2]),
            'posterior_probability_positive':float((beta>0).mean()),'max_rhat':rhat,
            'minimum_bulk_ess':ess,'divergences':div,'draws_per_chain':4000,'tune_per_chain':3000,
            'target_accept':.995,'status':'diagnostics_pass' if rhat<1.01 and ess>=400 and div==0 else 'WITHHOLD_DIAGNOSTICS_FAILED',
            'tree_sha256':sha256(tree),**design_info}
        results.append(result);(a.out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
        print(json.dumps(result),flush=True)
    (a.out/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')


if __name__=='__main__':main()
