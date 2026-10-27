"""Negative control: within the SAME habitat (animal), do source cohorts differ?
If Chicken_ vs Animal_ genomes of the same species differ in AMP count, then
cohort/batch effects of comparable size exist with habitat held constant."""
import sys, warnings
import numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
from scipy.stats import t as tdist
rng=np.random.default_rng(23)

g=pd.read_csv(sys.argv[1])
def resolved(l):
    n=l.fillna('').str.strip().str.removeprefix('s__').str.strip()
    return ~n.str.lower().isin(['','unknown','unclassified','uncultured'])
g=g[resolved(g.species)].copy()
g['cohort']=g.genome_id.str.replace(r'[_\-]?\d+$','',regex=True).str.replace(r'\d+','',regex=True)

print("=== aliasing check: is habitat a deterministic function of cohort? ===")
chk=g.groupby('cohort').habitat.nunique()
print(chk.to_string())
print(f"\ncohorts spanning >1 habitat: {(chk>1).sum()} of {len(chk)}")
print(f"habitats spanning >1 cohort: {g.groupby('habitat').cohort.nunique().to_dict()}")

BASE=['log_genome_bp','completeness','contamination','gc_fraction']
def fitpair(d,indicator,label,extra=()):
    d=d.copy()
    d['log_genome_bp']=np.log(d.genome_bp); d['log_contig_n50']=np.log(d.contig_n50)
    for c in BASE+[e for e in extra if e=='log_contig_n50']:
        s=d[c].std(ddof=0); d[c]=(d[c]-d[c].mean())/(s if s>0 else 1)
    terms=[indicator,'C(species)']+BASE+list(extra)
    f='amp_count ~ '+' + '.join(terms)
    m=smf.glm(f,data=d,family=sm.families.Poisson())
    if np.linalg.matrix_rank(m.exog)!=m.exog.shape[1]:
        print(f"{label:<46s} rank deficient"); return None
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r=m.fit(maxiter=300,cov_type='cluster',cov_kwds={'groups':d.species,'use_correction':True},use_t=True)
    b=float(r.params[indicator]); se=float(r.bse[indicator]); df=d.species.nunique()-1
    lo,hi=np.exp([b-tdist.ppf(.975,df)*se,b+tdist.ppf(.975,df)*se])
    print(f"{label:<46s} n={len(d):4d} k={d.species.nunique():2d}  "
          f"RR={np.exp(b):.3f} ({lo:.3f}-{hi:.3f})  p={2*tdist.sf(abs(b/se),df):.4g}")
    return np.exp(b)

print("\n=== NEGATIVE CONTROL: Chicken_ vs Animal_ , both habitat='animal' ===")
a=g[g.habitat.eq('animal') & g.cohort.isin(['Chicken','Animal'])].copy()
for high,tag in [(False,'published quality'),(True,'high quality (>=90/<5)')]:
    d=a[a.completeness.ge(90)&a.contamination.lt(5)].copy() if high else a.copy()
    d=d[d.groupby('species').cohort.transform('nunique').eq(2)].copy()
    d=d[d.groupby('species').amp_count.transform('sum').gt(0)].copy()
    if d.species.nunique()<3:
        print(f"  {tag}: only {d.species.nunique()} species span both cohorts -> cannot fit"); continue
    d['chicken']=d.cohort.eq('Chicken').astype(int)
    fitpair(d,'chicken',f"  Chicken/Animal within animal [{tag}]")
    fitpair(d,'chicken',f"  + log N50 [{tag}]",extra=('log_contig_n50',))
    print(f"     species spanning both cohorts: {d.species.nunique()}  "
          f"(chicken n={int(d.chicken.sum())}, animal n={int((1-d.chicken).sum())})")
    print(f"     zero-hit %: Chicken {100*(d[d.chicken==1].amp_count==0).mean():.1f}  "
          f"Animal {100*(d[d.chicken==0].amp_count==0).mean():.1f}")

print("\n=== raw cohort-level AMP rates (all resolved species) ===")
print(g.groupby(['habitat','cohort']).agg(genomes=('amp_count','size'),
    mean_amp=('amp_count','mean'),pct_zero=('amp_count',lambda s:100*(s==0).mean()),
    med_n50=('contig_n50','median'),med_compl=('completeness','median'),
    med_contam=('contamination','median'),med_bp=('genome_bp','median')).to_string())

print("\n=== positive-genome redundancy (clonality proxy, positives only) ===")
h=pd.read_csv(sys.argv[2])
pos=g[g.amp_count.gt(0)]
hh=h[h.genome_id.isin(pos.genome_id)]
prof=hh.groupby('genome_id').published_family.apply(lambda s:tuple(sorted(set(s))))
pp=pos.merge(prof.rename('profile'),left_on='genome_id',right_index=True,how='inner')
grp=pp.groupby(['species','habitat','profile']).size()
print(f"positive genomes: {len(pp)}   distinct species x habitat x family-profile: {len(grp)}")
print(f"redundancy among positives: {len(pp)/len(grp):.2f}x   groups >1: {(grp>1).sum()}  max {grp.max()}")
