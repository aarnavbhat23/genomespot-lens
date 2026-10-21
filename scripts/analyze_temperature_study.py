"""Fit prespecified REAL phylogenetic temperature/AMP count models.

No simulated fallback. Diagnostic failures withhold an interpretable primary
result, but fits/diagnostics remain available for repair and review.
"""
import argparse
from datetime import datetime, timezone
from importlib.metadata import version
import json
from pathlib import Path

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from genomespot_lens.amp_study import build_phylogenetic_nb, tree_covariance
from genomespot_lens.core import sha256

AA = "ACDEFGHIKLMNPQRSTVWY"


def design(data, temperature, adjusted):
    data = data.copy()
    counts = data[[f"aa_{a}" for a in AA]].to_numpy() * (data.reviewed_aa_fraction*data.proteome_aa).to_numpy()[:,None]
    logs = np.log(counts + .5);clr=logs-logs.mean(axis=1,keepdims=True)
    pca=PCA(n_components=3,svd_solver="full").fit(clr)
    x=pd.DataFrame({"Intercept":np.ones(len(data)),
        "temperature_per10C":(data[temperature].to_numpy()-data[temperature].mean())/10})
    values=np.c_[np.log(data.proteome_aa),np.log(data.reviewed_aa_fraction),pca.transform(clr)]
    scaled=StandardScaler().fit_transform(values)
    for i,name in enumerate(["log_proteome_aa_z","log_reviewed_fraction_z","aa_PC1_z","aa_PC2_z","aa_PC3_z"]):
        if i<2 or adjusted:x[name]=scaled[:,i]
    return x,{"aa_pca_variance_fraction":pca.explained_variance_ratio_.tolist(),
        "design_condition_number":float(np.linalg.cond(x)),"temperature_mean":float(data[temperature].mean()),
        "temperature_range":[float(data[temperature].min()),float(data[temperature].max())]}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    p.add_argument("--models",default="baseline,adjusted,parents100,broad,predicted")
    p.add_argument("--draws",type=int,default=1000);p.add_argument("--tune",type=int,default=1000)
    p.add_argument("--target-accept",type=float,default=.95)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    provenance=json.loads((a.data/"provenance.json").read_text())
    if provenance.get("data_kind")!="real":raise ValueError("Real input required")
    g=pd.read_csv(a.data/"genomes.csv")
    primary=g[g.primary_eligible].copy()
    if len(primary)<20:raise ValueError("Protocol feasibility threshold not met")
    models={"baseline":(primary,False,"observed_temperature"),
        "adjusted":(primary,True,"observed_temperature"),
        "parents100":(primary[primary.matched_source_proteins>=100].copy(),True,"observed_temperature"),
        "broad":(g[g.genome_qc_pass & g.observed_temperature.notna() & g.gtdb_tip.notna() & g.matched_source_proteins.gt(0)].copy(),True,"observed_temperature"),
        "predicted":(primary[primary.predicted_temperature_unflagged].copy(),True,"predicted_temperature")}
    results=[]
    for label in a.models.split(","):
        if label not in models:raise ValueError("Unknown model")
        data,adjusted,temperature=models[label];data=data.reset_index(drop=True)
        if len(data)<20:
            results.append({"model":label,"genomes":len(data),"status":"below_20_genome_threshold"});continue
        x,diagnostic=design(data,temperature,adjusted)
        cov=tree_covariance(a.data/"phylogeny.nwk",data.genome_id.tolist())
        model=build_phylogenetic_nb(data,x,cov)
        print('Fitting',label,'REAL genomes',len(data),flush=True)
        with model:
            fit=pm.sample(draws=a.draws,tune=a.tune,chains=4,cores=2,target_accept=a.target_accept,
                random_seed=23,progressbar=False,idata_kwargs={"log_likelihood":True})
            predictive=pm.sample_posterior_predictive(fit,random_seed=23,progressbar=False)
        fit.to_netcdf(a.out/f"{label}.nc")
        summary=az.summary(fit,var_names=["beta","phylo_sd","nb_shape"])
        summary.to_csv(a.out/f"{label}_diagnostics.csv")
        max_rhat=float(az.rhat(fit).to_array().max());min_ess=float(az.ess(fit,method="bulk").to_array().min())
        divergences=int(fit.sample_stats.diverging.sum())
        passed=max_rhat<1.01 and min_ess>=400 and divergences==0
        draws=fit.posterior.beta.sel(coefficient="temperature_per10C").values.flatten()
        rate=np.exp(draws);ci=np.quantile(rate,[.025,.5,.975])
        observed=data.amp_count.to_numpy();sim=predictive.posterior_predictive.amp_count.values.reshape(-1,len(data))
        checks={"observed_zero_count":int((observed==0).sum()),
            "posterior_predictive_zero_count_q025_q975":np.quantile((sim==0).sum(axis=1),[.025,.975]).tolist(),
            "observed_total":int(observed.sum()),"posterior_predictive_total_q025_q975":np.quantile(sim.sum(axis=1),[.025,.975]).tolist()}
        pd.DataFrame({"genome_id":data.genome_id,"observed_amp_count":observed,
            "predicted_median":np.median(sim,axis=0),"predicted_q025":np.quantile(sim,.025,axis=0),
            "predicted_q975":np.quantile(sim,.975,axis=0)}).to_csv(a.out/f"{label}_predictive.csv",index=False)
        result={"model":label,"genomes":len(data),"status":"diagnostics_pass" if passed else "WITHHOLD_DIAGNOSTICS_FAILED",
            "temperature_variable":temperature,"rr_per10C_median":float(ci[1]),"rr_q025":float(ci[0]),"rr_q975":float(ci[2]),
            "posterior_probability_positive":float((draws>0).mean()),"max_rhat":max_rhat,"minimum_bulk_ess":min_ess,
            "divergences":divergences,"draws_per_chain":a.draws,"tune_per_chain":a.tune,"target_accept":a.target_accept,
            **diagnostic,**checks}
        results.append(result)
        (a.out/f"{label}_summary.json").write_text(json.dumps(result,indent=2)+'\n')
        (a.out/"results.json").write_text(json.dumps(results,indent=2)+'\n')
        print(json.dumps(result),flush=True)
    (a.out/"run_provenance.json").write_text(json.dumps({"created_utc":datetime.now(timezone.utc).isoformat(),
        "data_kind":"real","data_provenance_sha256":sha256(a.data/"provenance.json"),
        "analysis_code_sha256":sha256(Path(__file__)),"models_requested":a.models,
        "packages":{n:version(n) for n in ['numpy','pandas','pymc','arviz','scikit-learn']}},indent=2)+'\n')


if __name__=="__main__":main()
