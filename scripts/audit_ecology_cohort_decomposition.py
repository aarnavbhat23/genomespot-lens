"""Decompose the human/animal lead by which animal source cohort supplies the comparator."""
import sys, warnings
import numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
from scipy.stats import t as tdist

g=pd.read_csv(sys.argv[1])
def resolved(l):
    n=l.fillna('').str.strip().str.removeprefix('s__').str.strip()
    return ~n.str.lower().isin(['','unknown','unclassified','uncultured'])
g=g[resolved(g.species)].copy()
g['cohort']=g.genome_id.str.replace(r'[_\-]?\d+$','',regex=True).str.replace(r'\d+','',regex=True)
BASE=['log_genome_bp','completeness','contamination','gc_fraction']

def fit(d,label,extra=()):
    d=d.copy()
    d['log_genome_bp']=np.log(d.genome_bp); d['log_contig_n50']=np.log(d.contig_n50)
    for c in BASE+[e for e in extra if e=='log_contig_n50']:
        s=d[c].std(ddof=0); d[c]=(d[c]-d[c].mean())/(s if s>0 else 1)
    terms=['habitat_indicator','C(species)']+BASE+list(extra)
    m=smf.glm('amp_count ~ '+' + '.join(terms),data=d,family=sm.families.Poisson())
    if np.linalg.matrix_rank(m.exog)!=m.exog.shape[1]:
        print(f"  {label:<44s} rank deficient"); return
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r=m.fit(maxiter=300,cov_type='cluster',cov_kwds={'groups':d.species,'use_correction':True},use_t=True)
    b=float(r.params['habitat_indicator']); se=float(r.bse['habitat_indicator'])
    k=d.species.nunique(); df=k-1
    lo,hi=np.exp([b-tdist.ppf(.975,df)*se,b+tdist.ppf(.975,df)*se])
    flag='' if k>=10 else '   [<10 clusters]'
    print(f"  {label:<44s} n={len(d):4d} k={k:2d}  RR={np.exp(b):.3f} ({lo:.3f}-{hi:.3f})  "
          f"p={2*tdist.sf(abs(b/se),df):.4g}{flag}")

def build(num_cohorts,den_cohorts,high=True):
    d=g[g.cohort.isin(num_cohorts+den_cohorts)].copy()
    d['habitat_indicator']=d.cohort.isin(num_cohorts).astype(int)
    if high: d=d[d.completeness.ge(90)&d.contamination.lt(5)].copy()
    d=d[d.groupby('species').habitat_indicator.transform('nunique').eq(2)].copy()
    d=d[d.groupby('species').amp_count.transform('sum').gt(0)].copy()
    return d

print("=== human/animal lead, decomposed by animal source cohort (high quality) ===")
print("  published frozen result: RR=0.555 (0.398-0.774), p=0.0018\n")
fit(build(['Human'],['Chicken','Animal']), "Human vs (Chicken + Animal)  [= frozen]")
fit(build(['Human'],['Chicken','Animal']), "Human vs (Chicken + Animal) + N50", ('log_contig_n50',))
print()
fit(build(['Human'],['Chicken']),          "Human vs Chicken only")
fit(build(['Human'],['Animal']),           "Human vs Animal only")
print()
print("=== same, published quality (no >=90/<5 filter) ===")
fit(build(['Human'],['Chicken','Animal'],high=False), "Human vs (Chicken + Animal)")
fit(build(['Human'],['Chicken'],high=False),          "Human vs Chicken only")
fit(build(['Human'],['Animal'],high=False),           "Human vs Animal only")

print("\n=== species overlap available for each comparison (high quality) ===")
for nm,nc,dc in [("Human vs Chicken",['Human'],['Chicken']),
                 ("Human vs Animal",['Human'],['Animal']),
                 ("Chicken vs Animal",['Chicken'],['Animal'])]:
    d=build(nc,dc)
    print(f"  {nm:<20s} species={d.species.nunique():2d}  genomes={len(d):4d}  "
          f"num={int(d.habitat_indicator.sum()):4d} den={int((1-d.habitat_indicator).sum()):4d}")

print("\n=== mean AMP count by cohort within shared species (high quality) ===")
d=g[g.completeness.ge(90)&g.contamination.lt(5)]
shared=d.groupby('species').cohort.nunique().pipe(lambda s:s[s>=2]).index
t=(d[d.species.isin(shared)].groupby(['species','cohort'])
     .agg(n=('amp_count','size'),mean=('amp_count','mean'),
          pct_zero=('amp_count',lambda s:100*(s==0).mean())).round(2))
print(t.to_string())
