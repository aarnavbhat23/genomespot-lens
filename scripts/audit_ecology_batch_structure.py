"""Reproduce the frozen human/animal lead, then add the controls and inference
the robustness protocol listed as untested: assembly fragmentation, flexible
genome-size adjustment, few-cluster-valid uncertainty."""
import sys, warnings
import numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
from scipy.stats import t as tdist

rng = np.random.default_rng(23)
g = pd.read_csv(sys.argv[1])

def resolved(labels):
    n = labels.fillna('').str.strip().str.removeprefix('s__').str.strip()
    return ~n.str.lower().isin(['', 'unknown', 'unclassified', 'uncultured'])

def cohort(reg, num, den, high=True):
    d = reg[reg.habitat.isin([num, den]) & resolved(reg.species)].copy()
    if high: d = d[d.completeness.ge(90) & d.contamination.lt(5)].copy()
    d = d[d.groupby('species').habitat.transform('nunique').eq(2)].copy()
    d = d[d.groupby('species').amp_count.transform('sum').gt(0)].copy().reset_index(drop=True)
    d['outcome'] = d.amp_count
    d['habitat_indicator'] = d.habitat.eq(num).astype(int)
    return d

BASE = ['log_genome_bp', 'completeness', 'contamination', 'gc_fraction']

def prep(d, extra=()):
    d = d.copy()
    d['log_genome_bp'] = np.log(d.genome_bp)
    d['log_contig_n50'] = np.log(d.contig_n50)
    d['log_genome_bp_sq'] = d['log_genome_bp'] ** 2
    for c in list(BASE) + [e for e in extra if e in
                           ('log_contig_n50', 'log_genome_bp_sq')]:
        s = d[c].std(ddof=0)
        d[c] = (d[c] - d[c].mean()) / (s if s > 0 else 1.0)
    return d

def fit(d, extra=(), family=None):
    terms = ['habitat_indicator', 'C(species)'] + BASE + list(extra)
    f = 'outcome ~ ' + ' + '.join(terms)
    fam = family or sm.families.Poisson()
    m = smf.glm(f, data=d, family=fam)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r = m.fit(maxiter=300, cov_type='cluster',
                  cov_kwds={'groups': d.species, 'use_correction': True}, use_t=True)
    b = float(r.params['habitat_indicator']); se = float(r.bse['habitat_indicator'])
    df = d.species.nunique() - 1
    lo, hi = np.exp([b - tdist.ppf(.975, df) * se, b + tdist.ppf(.975, df) * se])
    return dict(rr=np.exp(b), lo=lo, hi=hi,
                p=float(2 * tdist.sf(abs(b / se), df)), beta=b, se=se,
                n=len(d), k=d.species.nunique(),
                disp=float(r.pearson_chi2 / r.df_resid)), r, m

def wild_cluster_bootstrap(d, extra=(), reps=9999):
    """Rademacher wild cluster bootstrap-t, imposing the null (Cameron-Gelbach-Miller).
    Appropriate when the number of clusters is small (here 16)."""
    terms = ['C(species)'] + BASE + list(extra)
    f_null = 'outcome ~ ' + ' + '.join(terms)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        null = smf.glm(f_null, data=d, family=sm.families.Poisson()).fit(maxiter=300)
    mu0 = np.asarray(null.fittedvalues)
    resid = d.outcome.values - mu0
    obs, _, _ = fit(d, extra)
    tobs = abs(obs['beta'] / obs['se'])
    sp = d.species.values; uniq = np.unique(sp)
    cnt = 0; ok = 0
    for _ in range(reps):
        w = pd.Series(rng.choice([-1.0, 1.0], size=len(uniq)), index=uniq).reindex(sp).values
        y = np.round(np.clip(mu0 + w * resid, 0, None)).astype(int)
        db = d.copy(); db['outcome'] = y
        if db.groupby('species').outcome.sum().gt(0).sum() < db.species.nunique():
            continue
        try:
            r, _, _ = fit(db, extra)
            if not np.isfinite(r['se']) or r['se'] <= 0: continue
            ok += 1
            if abs(r['beta'] / r['se']) >= tobs: cnt += 1
        except Exception:
            continue
    return (cnt + 1) / (ok + 1), ok

def randomization_within_species(d, extra=(), reps=9999):
    """Permute habitat labels within species: exact-style test of the
    within-species estimand, free of asymptotic cluster assumptions."""
    obs, _, _ = fit(d, extra)
    tobs = abs(obs['beta'] / obs['se'])
    cnt = 0; ok = 0
    for _ in range(reps):
        db = d.copy()
        db['habitat_indicator'] = (db.groupby('species')['habitat_indicator']
                                   .transform(lambda s: rng.permutation(s.values)))
        if db.habitat_indicator.nunique() < 2: continue
        try:
            r, _, _ = fit(db, extra)
            if not np.isfinite(r['se']) or r['se'] <= 0: continue
            ok += 1
            if abs(r['beta'] / r['se']) >= tobs: cnt += 1
        except Exception:
            continue
    return (cnt + 1) / (ok + 1), ok

def show(tag, r):
    print(f"{tag:<52s} n={r['n']:4d} k={r['k']:2d}  RR={r['rr']:.3f} "
          f"({r['lo']:.3f}-{r['hi']:.3f})  p={r['p']:.4g}  disp={r['disp']:.2f}")

d = prep(cohort(g, 'human', 'animal', high=True),
         extra=('log_contig_n50', 'log_genome_bp_sq'))
print(f"=== high-quality human/animal count cohort: {len(d)} genomes / "
      f"{d.species.nunique()} species ===\n")

print("--- A. reproduce frozen specification ---")
base, bfit, _ = fit(d); show("frozen: species + size + compl + contam + GC", base)
print("   (published: RR=0.555, 0.398-0.774, p=0.00182)\n")

print("--- B. add untested assembly-fragmentation controls ---")
for tag, extra in [("+ log contig N50", ('log_contig_n50',)),
                   ("+ log N50 + quadratic genome size", ('log_contig_n50', 'log_genome_bp_sq')),
                   ("+ quadratic genome size only", ('log_genome_bp_sq',)),
                   ("+ assembly status", ('C(assembly_status)',))]:
    try:
        r, _, _ = fit(d, extra); show(tag, r)
    except Exception as e:
        print(f"{tag:<52s} FAILED {type(e).__name__}: {e}")

print("\n--- C. negative binomial (dispersion check) ---")
for a in (0.1, 0.5, 1.0):
    try:
        r, _, _ = fit(d, (), sm.families.NegativeBinomial(alpha=a))
        show(f"NB2 alpha={a}", r)
    except Exception as e:
        print(f"NB2 alpha={a} FAILED {e}")

print("\n--- D. few-cluster-valid inference on the frozen specification ---")
print("    (16 species clusters; sandwich t intervals are known to under-cover)")
p_wcb, n_wcb = wild_cluster_bootstrap(d, (), reps=1999)
print(f"    wild cluster bootstrap-t (Rademacher, null imposed): p={p_wcb:.4f}  [{n_wcb} usable reps]")
p_ri, n_ri = randomization_within_species(d, (), reps=1999)
print(f"    within-species randomization test:                   p={p_ri:.4f}  [{n_ri} usable reps]")

print("\n--- E. fragmentation / size differ by habitat? (confounder check) ---")
for c in ['contig_n50', 'genome_bp', 'completeness', 'contamination', 'gc_fraction']:
    raw = cohort(g, 'human', 'animal', high=True)
    h = raw[raw.habitat.eq('human')][c]; a = raw[raw.habitat.eq('animal')][c]
    # within-species demeaned difference, matching the estimand
    dm = raw.assign(v=raw[c]).copy()
    dm['v_c'] = dm.v - dm.groupby('species').v.transform('mean')
    hd = dm[dm.habitat.eq('human')].v_c.mean(); ad = dm[dm.habitat.eq('animal')].v_c.mean()
    print(f"  {c:<14s} human med={h.median():>12.4g}  animal med={a.median():>12.4g}   "
          f"within-species mean diff (h-a)={hd-ad:>12.4g}")
