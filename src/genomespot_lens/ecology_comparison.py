"""Species-controlled count comparisons for a fully enumerated genome registry."""
import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.api as sm
import statsmodels.formula.api as smf

HABITATS = ['animal', 'human', 'fermented food']
TERMS = {h: f"C(habitat, Treatment(reference='animal'))[T.{h}]" for h in HABITATS[1:]}


def fit_count_contrasts(data, *, species_effect=False, quality=False, label='crude'):
    """PPML with species-cluster sandwich uncertainty; counts are not potency."""
    d = data.copy()
    dropped_zero = 0
    if species_effect:
        positive = d.groupby('species').amp_count.transform('sum').gt(0)
        dropped_zero = int((~positive).sum())
        d = d[positive].copy()
    formula = "amp_count ~ C(habitat, Treatment(reference='animal'))"
    if species_effect:
        formula += ' + C(species)'
    if quality:
        d['log_genome_bp'] = np.log(d.genome_bp)
        d['completeness_per10'] = d.completeness/10
        for name in ['log_genome_bp','completeness_per10','contamination','gc_fraction']:
            scale = d[name].std(ddof=0)
            if not scale > 0: raise ValueError('Constant quality covariate')
            d[name] = (d[name]-d[name].mean())/scale
        formula += ' + log_genome_bp + completeness_per10 + contamination + gc_fraction'
    model = smf.glm(formula, data=d, family=sm.families.Poisson())
    if np.linalg.matrix_rank(model.exog) != model.exog.shape[1]:
        raise ValueError('Rank deficient comparison')
    groups = d.species.nunique()
    if groups < 10:
        raise ValueError('Too few independent species clusters')
    fit = model.fit(maxiter=200, cov_type='cluster',
        cov_kwds={'groups':d.species,'use_correction':True}, use_t=True)
    if not fit.converged:
        raise ValueError('Count model did not converge')
    names = fit.params.index.tolist(); covariance = fit.cov_params().to_numpy()
    habitat_indices = [names.index(TERMS[h]) for h in HABITATS[1:]]
    b = fit.params.to_numpy()[habitat_indices]
    v = covariance[np.ix_(habitat_indices,habitat_indices)]
    omnibus_f = float(b @ np.linalg.solve(v,b)/2)
    omnibus_p = float(st.f.sf(omnibus_f,2,groups-1))
    rows = []
    for numerator, denominator in [('human','animal'),('fermented food','animal'),('fermented food','human')]:
        contrast = np.zeros(len(names))
        if numerator != 'animal': contrast[names.index(TERMS[numerator])] += 1
        if denominator != 'animal': contrast[names.index(TERMS[denominator])] -= 1
        beta = float(contrast @ fit.params.to_numpy())
        variance = float(contrast @ covariance @ contrast)
        if variance <= 0: raise ValueError('Nonpositive contrast variance')
        se = np.sqrt(variance); df = groups-1; critical = st.t.ppf(.975,df)
        lo,hi = np.exp([beta-critical*se, beta+critical*se])
        rows.append({'model':label,'contrast':f'{numerator} / {denominator}',
            'rate_ratio':float(np.exp(beta)),'ci_low':float(lo),'ci_high':float(hi),
            'p_value':float(2*st.t.sf(abs(beta/se),df)),
            'omnibus_habitat_f':omnibus_f,'omnibus_habitat_p':omnibus_p,
            'genomes':len(d),'species_clusters':groups,'df':df,
            'all_zero_species_genomes_conditioned_out':dropped_zero,
            'species_fixed_effect':species_effect,'quality_adjusted':quality,
            'converged':bool(fit.converged),
            'pearson_dispersion':float(fit.pearson_chi2/fit.df_resid),
            'design_condition_number':float(np.linalg.cond(model.exog))})
    return pd.DataFrame(rows)


def within_species_permutation(data, permutations=9999, seed=23):
    """Omnibus species-conditioned Poisson score, permuting habitat within species.

    No quality adjustment. Exchangeability within species is an assumption;
    source-study, fine-scale strain structure and quality can still confound.
    The information matrix is unchanged by label permutations because species
    totals and within-species habitat counts are held fixed.
    """
    d = data.reset_index(drop=True)
    species, _ = pd.factorize(d.species)
    habitat = pd.Categorical(d.habitat,categories=HABITATS).codes
    if (habitat < 0).any(): raise ValueError('Unknown habitat')
    y = d.amp_count.to_numpy(float)
    sizes = np.bincount(species)
    totals = np.bincount(species,weights=y)
    expected = totals/sizes
    residual = y-expected[species]
    x = np.c_[habitat==1,habitat==2].astype(float)
    centered = x.copy()
    for j in range(2):
        centered[:,j] -= (np.bincount(species,weights=x[:,j])/sizes)[species]
    information = centered.T @ (centered*expected[species,None])
    if np.linalg.matrix_rank(information) != 2: raise ValueError('No full habitat contrast')
    inverse = np.linalg.inv(information)
    def statistic(labels):
        score = np.array([residual[labels==1].sum(),residual[labels==2].sum()])
        return float(score @ inverse @ score)
    observed = statistic(habitat)
    indices = [np.flatnonzero(species==s) for s in range(len(sizes))]
    varying = [ix for ix in indices if len(set(habitat[ix]))>1]
    rng = np.random.default_rng(seed); exceed = 0
    for _ in range(permutations):
        shuffled = habitat.copy()
        for ix in varying: shuffled[ix] = rng.permutation(habitat[ix])
        exceed += statistic(shuffled) >= observed-1e-12
    return {'statistic':observed,'permutations':permutations,'seed':seed,
        'exceedances':exceed,'p_value':(exceed+1)/(permutations+1),
        'genomes':len(d),'species':len(sizes),'multi_habitat_species':len(varying),
        'quality_adjusted':False,
        'null_assumption':'Habitat labels exchangeable within species; observational sensitivity only'}
