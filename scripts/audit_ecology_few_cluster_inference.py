import sys, warnings
import numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
from scipy.stats import t as tdist
rng=np.random.default_rng(23)
g=pd.read_csv(sys.argv[1])
def resolved(l):
    n=l.fillna('').str.strip().str.removeprefix('s__').str.strip()
    return ~n.str.lower().isin(['','unknown','unclassified','uncultured'])
BASE=['log_genome_bp','completeness','contamination','gc_fraction']
def cohort(reg,num,den,high=True):
    d=reg[reg.habitat.isin([num,den]) & resolved(reg.species)].copy()
    if high: d=d[d.completeness.ge(90)&d.contamination.lt(5)].copy()
    d=d[d.groupby('species').habitat.transform('nunique').eq(2)].copy()
    d=d[d.groupby('species').amp_count.transform('sum').gt(0)].copy().reset_index(drop=True)
    d['habitat_indicator']=d.habitat.eq(num).astype(int); return d
def prep(d):
    d=d.copy(); d['log_genome_bp']=np.log(d.genome_bp)
    for c in BASE:
        s=d[c].std(ddof=0); d[c]=(d[c]-d[c].mean())/(s if s>0 else 1)
    return d
def fitb(d):
    m=smf.glm('amp_count ~ habitat_indicator + C(species) + '+' + '.join(BASE),
              data=d,family=sm.families.Poisson())
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r=m.fit(maxiter=300,cov_type='cluster',
                cov_kwds={'groups':d.species,'use_correction':True},use_t=True)
    return float(r.params['habitat_indicator']),float(r.bse['habitat_indicator'])

d=prep(cohort(g,'human','animal'))
print("=== 1. design imbalance in the 16 informative species ===")
tab=d.groupby('species').agg(human=('habitat_indicator','sum'),total=('habitat_indicator','size'))
tab['animal']=tab.total-tab.human
tab['min_arm']=tab[['human','animal']].min(axis=1)
tab['amp_total']=d.groupby('species').amp_count.sum()
print(tab.sort_values('min_arm').to_string())
print(f"\nspecies with an arm of size 1 : {(tab.min_arm==1).sum()}")
print(f"species with an arm <= 3      : {(tab.min_arm<=3).sum()}")
print(f"sum of min-arm sizes (rough effective paired information): {tab.min_arm.sum()} of {len(d)} genomes")
print(f"total AMP hits in whole cohort: {int(d.amp_count.sum())}")

print("\n=== 2. stable few-cluster inference, frozen specification (9999 reps) ===")
b,se=fitb(d); tobs=abs(b/se)
print(f"    point: RR={np.exp(b):.3f}  sandwich p={2*tdist.sf(tobs,d.species.nunique()-1):.5f}")
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    null=smf.glm('amp_count ~ C(species) + '+' + '.join(BASE),data=d,
                 family=sm.families.Poisson()).fit(maxiter=300)
mu0=np.asarray(null.fittedvalues); resid=d.amp_count.values-mu0
sp=d.species.values; uniq=np.unique(sp)
cnt=ok=0
for _ in range(9999):
    w=pd.Series(rng.choice([-1.,1.],size=len(uniq)),index=uniq).reindex(sp).values
    y=np.round(np.clip(mu0+w*resid,0,None)).astype(int)
    db=d.copy(); db['amp_count']=y
    if db.groupby('species').amp_count.sum().gt(0).sum()<db.species.nunique(): continue
    try:
        bb,ss=fitb(db)
        if not np.isfinite(ss) or ss<=0: continue
        ok+=1
        if abs(bb/ss)>=tobs: cnt+=1
    except Exception: continue
print(f"    wild cluster bootstrap-t p = {(cnt+1)/(ok+1):.4f}   [{ok} usable reps]")
cnt=ok=0
for _ in range(9999):
    db=d.copy()
    db['habitat_indicator']=db.groupby('species')['habitat_indicator'].transform(
        lambda s: rng.permutation(s.values))
    if db.habitat_indicator.nunique()<2: continue
    try:
        bb,ss=fitb(db)
        if not np.isfinite(ss) or ss<=0: continue
        ok+=1
        if abs(bb/ss)>=tobs: cnt+=1
    except Exception: continue
print(f"    within-species randomization p = {(cnt+1)/(ok+1):.4f}   [{ok} usable reps]")

print("\n=== 3. leave-one-species-out on the frozen specification ===")
out=[]
for s in sorted(d.species.unique()):
    dd=prep(cohort(g,'human','animal'))
    dd=dd[dd.species!=s].copy()
    dd=dd[dd.groupby('species').habitat_indicator.transform('nunique').eq(2)]
    dd=dd[dd.groupby('species').amp_count.transform('sum').gt(0)]
    try:
        bb,ss=fitb(dd); k=dd.species.nunique(); df=k-1
        lo,hi=np.exp([bb-tdist.ppf(.975,df)*ss,bb+tdist.ppf(.975,df)*ss])
        out.append({'dropped':s.removeprefix('s__'),'k':k,'RR':np.exp(bb),
                    'lo':lo,'hi':hi,'crosses_1':bool(lo<=1<=hi)})
    except Exception as e:
        out.append({'dropped':s.removeprefix('s__'),'k':np.nan,'RR':np.nan,
                    'lo':np.nan,'hi':np.nan,'crosses_1':None})
o=pd.DataFrame(out).sort_values('RR')
print(o.round(3).to_string(index=False))
print(f"\nomissions whose interval includes 1: {o.crosses_1.sum()} of {len(o)}")
