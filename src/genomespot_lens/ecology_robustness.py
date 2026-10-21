"""Direct habitat contrasts with explicit common support and endpoint separation."""
import warnings
import numpy as np
import pandas as pd
from scipy.stats import t
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import PerfectSeparationError, PerfectSeparationWarning

PAIRS=[('human','animal'),('fermented food','animal'),('fermented food','human')]


def resolved_species(labels):
    names=labels.fillna('').str.strip().str.removeprefix('s__').str.strip()
    return ~names.str.lower().isin(['','unknown','unclassified','uncultured'])


def pair_cohort(registry,numerator,denominator,*,high_quality=False,presence=False):
    d=registry[registry.habitat.isin([numerator,denominator]) & resolved_species(registry.species)].copy()
    if high_quality:d=d[d.completeness.ge(90)&d.contamination.lt(5)].copy()
    overlap=d.groupby('species').habitat.transform('nunique').eq(2)
    d=d[overlap].copy()
    d['outcome']=(d.amp_count.gt(0).astype(int) if presence else d.amp_count)
    if presence:informative=d.groupby('species').outcome.transform('nunique').eq(2)
    else:informative=d.groupby('species').outcome.transform('sum').gt(0)
    d=d[informative].copy().reset_index(drop=True)
    d['habitat_indicator']=d.habitat.eq(numerator).astype(int)
    return d


def fit_direct_pair(registry,numerator,denominator,*,high_quality=False,presence=False):
    d=pair_cohort(registry,numerator,denominator,high_quality=high_quality,presence=presence)
    result={'contrast':f'{numerator} / {denominator}',
        'endpoint':'presence' if presence else 'count',
        'quality':'high' if high_quality else 'published',
        'effect_measure':'odds_ratio' if presence else 'count_rate_ratio',
        'genomes':len(d),'species_clusters':d.species.nunique(),
        'genome_ids':d.genome_id.tolist(),'species_ids':sorted(d.species.unique().tolist()),
        'status':'insufficient_support'}
    if result['species_clusters']<10:return result
    # Scaling changes numerical conditioning, not the specified mean model.
    d['log_genome_bp']=np.log(d.genome_bp)
    for name in ['log_genome_bp','completeness','contamination','gc_fraction']:
        scale=d[name].std(ddof=0)
        if not scale>0:
            result['status']='constant_covariate';return result
        d[name]=(d[name]-d[name].mean())/scale
    formula='outcome ~ habitat_indicator + C(species) + log_genome_bp + completeness + contamination + gc_fraction'
    model=smf.glm(formula,data=d,family=sm.families.Binomial() if presence else sm.families.Poisson())
    if np.linalg.matrix_rank(model.exog)!=model.exog.shape[1]:
        result['status']='rank_deficient';return result
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('error',PerfectSeparationWarning)
            fit=model.fit(maxiter=200,cov_type='cluster',cov_kwds={'groups':d.species,'use_correction':True},use_t=True)
    except (PerfectSeparationError,PerfectSeparationWarning,np.linalg.LinAlgError) as e:
        result['status']='withhold_separation_or_covariance';result['reason']=str(e);return result
    result.update({'converged':bool(fit.converged),'maximum_abs_coefficient':float(np.abs(fit.params).max())})
    if not fit.converged or not np.isfinite(fit.params).all() or np.abs(fit.params).max()>25:
        result['status']='withhold_nonconverged_or_extreme';return result
    beta=float(fit.params['habitat_indicator']);se=float(fit.bse['habitat_indicator'])
    if not np.isfinite(se) or se<=0:
        result['status']='withhold_invalid_variance';return result
    df=d.species.nunique()-1;lo,hi=np.exp([beta-t.ppf(.975,df)*se,beta+t.ppf(.975,df)*se])
    result.update({'status':'fit_pass','effect':float(np.exp(beta)),'ci_low':float(lo),'ci_high':float(hi),
        'p_value':float(2*t.sf(abs(beta/se),df)),'df':df,'converged':bool(fit.converged),
        'design_condition_number':float(np.linalg.cond(model.exog)),
        'pearson_dispersion':float(fit.pearson_chi2/fit.df_resid)})
    return result
