"""Exploratory prior and tree sensitivity; preserve the frozen original fits."""
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
from genomespot_lens.amp_study import tree_covariance
from genomespot_lens.core import sha256


def prune_exact_tree(source, destination, genome_ids):
    """Prune exact versioned tips, collapsing unary paths without changing distances."""
    with gzip.open(source, 'rt') as fh:
        tree = Phylo.read(fh, 'newick')
    aliases = {t.name: t.name.removeprefix('RS_').removeprefix('GB_')
               for t in tree.get_terminals()}
    selected = {tip: asm for tip, asm in aliases.items() if asm in set(genome_ids)}
    if len(set(selected.values())) != len(selected):
        raise ValueError('Ambiguous exact assembly tip')
    def retain(node):
        if node.is_terminal():
            if node.name not in selected:
                return None
            node.name = selected[node.name]
            return node
        node.clades = [r for c in node.clades if (r := retain(c)) is not None]
        if not node.clades:
            return None
        if len(node.clades) == 1:
            child = node.clades[0]
            child.branch_length = (child.branch_length or 0) + (node.branch_length or 0)
            return child
        return node
    tree.root = retain(tree.root)
    if tree.root is None:
        raise ValueError('No exact tree tips')
    tree.root.branch_length = 0
    tree.rooted = True
    Phylo.write(tree, destination, 'newick', format_branch_length='%1.10f')
    return set(selected.values())


def sensitivity_model(data, x, covariance, slope_sd):
    scales = np.full(x.shape[1], slope_sd)
    scales[list(x.columns).index('Intercept')] = 1.5
    if np.linalg.matrix_rank(x) != x.shape[1]:
        raise ValueError('Rank deficient design')
    chol = np.linalg.cholesky(covariance + np.eye(len(data)) * 1e-8)
    with pm.Model(coords={'genome': data.genome_id.tolist(), 'coefficient': list(x.columns)}) as model:
        beta = pm.Normal('beta', 0, scales, dims='coefficient')
        sd = pm.HalfNormal('phylo_sd', 1)
        z = pm.Normal('phylo_z', 0, 1, dims='genome')
        u = pm.Deterministic('phylo_effect', sd * pm.math.dot(chol, z), dims='genome')
        shape = pm.Exponential('nb_shape', 1)
        mu = pm.math.exp(np.log(data.n_screened.to_numpy() / 100000) + pm.math.dot(x, beta) + u)
        pm.NegativeBinomial('amp_count', mu=mu, alpha=shape,
                           observed=data.amp_count.to_numpy().astype(int), dims='genome')
    return model


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data', type=Path, required=True)
    p.add_argument('--tree-r232', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--models', default='prior075,prior300,tree232')
    p.add_argument('--draws', type=int, default=4000)
    p.add_argument('--tune', type=int, default=3000)
    a = p.parse_args(); a.out.mkdir(parents=True, exist_ok=True)
    if json.loads((a.data/'provenance.json').read_text())['data_kind'] != 'real':
        raise ValueError('Real input required')
    g = pd.read_csv(a.data/'genomes.csv')
    g = g[g.primary_eligible].reset_index(drop=True)
    common = prune_exact_tree(a.tree_r232, a.out/'phylogeny-r232.nwk', g.genome_id)
    specifications = {
        'prior075': (g, a.data/'phylogeny.nwk', .75),
        'prior300': (g, a.data/'phylogeny.nwk', 3.),
        'tree232': (g[g.genome_id.isin(common)].reset_index(drop=True), a.out/'phylogeny-r232.nwk', 1.5),
        'tree226_common': (g[g.genome_id.isin(common)].reset_index(drop=True), a.data/'phylogeny.nwk', 1.5),
    }
    requested = a.models.split(',')
    if common != set(g.genome_id) and 'tree232' in requested:
        requested.append('tree226_common')
    provenance = {'created_utc': datetime.now(timezone.utc).isoformat(),
        'protocol_commit': '07a490c', 'data_kind': 'real',
        'analysis_code_sha256': sha256(Path(__file__)),
        'original_design_code_sha256': sha256(Path(__file__).with_name('analyze_temperature_study.py')),
        'data_sha256': {n: sha256(a.data/n) for n in ['provenance.json','genomes.csv','phylogeny.nwk']},
        'r232_source_url': 'https://data.gtdb.ecogenomic.org/releases/release232/232.0/ar53_r232.tree.gz',
        'r232_source_sha256': sha256(a.tree_r232),
        'r232_pruned_sha256': sha256(a.out/'phylogeny-r232.nwk'),
        'primary_genome_ids': g.genome_id.tolist(), 'r232_common_genome_ids': sorted(common),
        'rooting': 'Published release root; exact versioned tips; unary paths collapsed',
        'pymc_version': pm.__version__, 'arviz_version': az.__version__}
    (a.out/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    results = []
    for label in requested:
        data, tree, prior = specifications[label]
        if len(data) < 20:
            raise ValueError('Below original 20-genome threshold')
        x, diagnostic = design(data, 'observed_temperature', True)
        cov = tree_covariance(tree, data.genome_id.tolist())
        model = sensitivity_model(data, x, cov, prior)
        print('Fitting', label, 'n=', len(data), flush=True)
        with model:
            fit = pm.sample(draws=a.draws, tune=a.tune, chains=4, cores=2,
                target_accept=.995, random_seed=23, progressbar=False,
                idata_kwargs={'log_likelihood': True})
        fit.to_netcdf(a.out/f'{label}.nc')
        az.summary(fit, var_names=['beta','phylo_sd','nb_shape']).to_csv(a.out/f'{label}_diagnostics.csv')
        rhat = float(az.rhat(fit).to_array().max())
        ess = float(az.ess(fit, method='bulk').to_array().min())
        div = int(fit.sample_stats.diverging.sum())
        beta = fit.posterior.beta.sel(coefficient='temperature_per10C').values.ravel()
        ci = np.quantile(np.exp(beta), [.025,.5,.975])
        result = {'model': label, 'genomes': len(data), 'non_intercept_prior_sd': prior,
            'rr_q025': float(ci[0]), 'rr_per10C_median': float(ci[1]), 'rr_q975': float(ci[2]),
            'posterior_probability_positive': float((beta>0).mean()),
            'max_rhat': rhat, 'minimum_bulk_ess': ess, 'divergences': div,
            'draws_per_chain': a.draws, 'tune_per_chain': a.tune, 'target_accept': .995,
            'status': 'diagnostics_pass' if rhat<1.01 and ess>=400 and div==0 else 'WITHHOLD_DIAGNOSTICS_FAILED',
            'genome_ids': data.genome_id.tolist(), **diagnostic}
        results.append(result)
        (a.out/'results.json').write_text(json.dumps(results, indent=2)+'\n')
        print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
