"""Real within-homology sensitivity with two-way cluster-robust uncertainty."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import patsy
from scipy.stats import t
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
from statsmodels.stats.sandwich_covariance import cov_cluster_2groups
from threadpoolctl import threadpool_limits

from genomespot_lens.core import sha256

AA="ACDEFGHIKLMNPQRSTVWY"


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data",type=Path,required=True);parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args();a.out.mkdir(parents=True,exist_ok=False)
    if json.loads((a.data/"provenance.json").read_text())["data_kind"]!="real":raise ValueError("Real data required")
    g=pd.read_csv(a.data/"genomes.csv");g=g[g.primary_eligible]
    proteins=pd.read_csv(a.data/"proteins.csv");proteins=proteins[proteins.genome_id.isin(g.genome_id)&proteins.n_screened.gt(0)].copy()
    proteins=proteins.merge(g[["genome_id","observed_temperature"]],on="genome_id",validate="m:1")
    logs=np.log(proteins[[f"aa_{aa}" for aa in AA]].to_numpy()*proteins.length_aa.to_numpy()[:,None]+.5)
    clr=logs-logs.mean(axis=1,keepdims=True);pca=PCA(n_components=3,svd_solver="full").fit(clr)
    pcs=StandardScaler().fit_transform(pca.transform(clr))
    for i in range(3):proteins[f'pc{i+1}']=pcs[:,i]
    proteins['temp10']=(proteins.observed_temperature-proteins.observed_temperature.mean())/10
    proteins['log_length_z']=StandardScaler().fit_transform(np.log(proteins[['length_aa']]))
    results=[]
    for cutoff in [30,50]:
        col=f'homology{cutoff}';coverage=proteins.groupby(col).genome_id.nunique();keep=coverage[coverage>=5].index
        d=proteins[proteins[col].isin(keep)].copy();totals=d.groupby(col).amp_count.sum();informative=totals[totals>0].index
        dropped=len(totals)-len(informative);d=d[d[col].isin(informative)].copy();d['family']=d[col]
        if len(informative)<10:raise ValueError('Fewer than ten informative homolog groups')
        for adjusted in [False,True]:
            formula='1 + temp10 + log_length_z + C(family)' + (' + pc1 + pc2 + pc3' if adjusted else '')
            x=patsy.dmatrix(formula,d,return_type='dataframe')
            if np.linalg.matrix_rank(x)!=x.shape[1]:raise ValueError('Rank-deficient within-homology design')
            with threadpool_limits(limits=1):
                fit=sm.GLM(d.amp_count,x,family=sm.families.Poisson(),offset=np.log(d.n_screened)).fit(maxiter=100)
                covariance,_,_=cov_cluster_2groups(fit,pd.Categorical(d.genome_id).codes,pd.Categorical(d.family).codes)
            i=x.columns.get_loc('temp10');variance=float(covariance[i,i]);beta=float(fit.params.iloc[i])
            df=min(d.genome_id.nunique(),d.family.nunique())-1
            se=np.sqrt(variance) if variance>0 else np.nan
            critical=t.ppf(.975,df)
            row={'identity_cutoff':cutoff,'composition_adjusted':adjusted,'parent_rows':len(d),'genomes':d.genome_id.nunique(),
                'informative_homology_groups':d.family.nunique(),'zero_total_groups_conditioned_out':dropped,
                'rr_per10C':float(np.exp(beta)),'rr_q025':float(np.exp(beta-critical*se)),
                'rr_q975':float(np.exp(beta+critical*se)),'pvalue':float(2*t.sf(abs(beta/se),df)),
                'cluster_degrees_freedom':df,'converged':bool(fit.converged),
                'temperature_variance_positive':variance>0,'condition_number':float(np.linalg.cond(x)),
                'pearson_dispersion':float(fit.pearson_chi2/fit.df_resid),
                'interpretation':'Secondary quasi-Poisson family contrast; two-way cluster robust, not phylogenetic mixed model'}
            results.append(row);print(json.dumps(row),flush=True)
    pd.DataFrame(results).to_csv(a.out/'effects.csv',index=False)
    (a.out/'provenance.json').write_text(json.dumps({'data_kind':'real','analysis_code_sha256':sha256(Path(__file__)),
        'data_sha256':sha256(a.data/'provenance.json'),'aa_pca_variance_fraction':pca.explained_variance_ratio_.tolist(),
        'sensitivity_only':True,'primary_genomes':len(g)},indent=2)+'\n')


if __name__=='__main__':main()
