import sys
import pandas as pd, numpy as np

df = pd.read_csv(sys.argv[1])
mic = pd.read_csv(sys.argv[2], index_col=0)

df['core'] = df.assembly.str.replace(r'^GC[AF]_', '', regex=True)
df['db']   = df.assembly.str[:3]

print("=== GCA/GCF redundancy ===")
per_core = df.groupby('core').db.nunique()
print(f"distinct accession strings      : {df.assembly.nunique()}")
print(f"distinct numeric cores          : {df.core.nunique()}")
print(f"cores present as BOTH GCA & GCF : {(per_core==2).sum()}")
print(f"cores GCA-only or GCF-only      : {(per_core==1).sum()}")

# Are the paired GCA/GCF entries identical in content?
both = per_core[per_core==2].index
chk = df[df.core.isin(both)].groupby(['core','db']).apply(
    lambda d: tuple(sorted(d.seq)), include_groups=False)
ident = chk.unstack('db').apply(lambda r: r['GCA'] == r['GCF'], axis=1)
print(f"paired cores with identical peptide sets: {ident.sum()} / {len(ident)}")

print("\n=== peptide redundancy ===")
occ = df.groupby('seq').size().sort_values(ascending=False)
print(f"header records            : {len(df)}")
print(f"distinct peptide sequences: {len(occ)}")
print(f"occurrences per peptide   : median {occ.median():.0f}, mean {occ.mean():.1f}, max {occ.max()}")
print(f"peptides found in 1 record only: {(occ==1).sum()}")
print("\nmost widespread peptides (records, distinct cores):")
top = occ.head(8).index
for s in top:
    sub=df[df.seq==s]
    print(f"  {sub.core.nunique():5d} cores  {len(sub):5d} records  {s[:45]}")

# dedup to core level
dedup = df.drop_duplicates(subset=['core','seq'])
print(f"\nrecords after GCA/GCF dedup (core+seq): {len(dedup)}")
print(f"per-core candidate counts: median {dedup.groupby('core').size().median():.0f}, "
      f"mean {dedup.groupby('core').size().mean():.2f}, max {dedup.groupby('core').size().max()}")

print("\n=== selection-rule sensitivity on released predicted MICs ===")
mn = mic.min(axis=1); md = mic.median(axis=1)
rules = {
    'any target  < 80 uM' : (mn <  80).sum(),
    'any target <= 80 uM' : (mn <= 80).sum(),
    'median     <= 80 uM' : (md <= 80).sum(),
    'median      < 80 uM' : (md <  80).sum(),
    'any target  < 40 uM' : (mn <  40).sum(),
    'all targets<= 80 uM' : (mic.max(axis=1) <= 80).sum(),
}
for k,v in rules.items():
    print(f"  {k:22s} {v:4d} / {len(mic)}  ({100*v/len(mic):5.1f}%)")
print(f"\n  ratio any<80 : median<=80  = {rules['any target  < 80 uM']/rules['median     <= 80 uM']:.2f}x")

print("\n=== MIC table coverage vs released FASTAs ===")
cand_seqs=set(df.seq); micseqs=set(mic.index)
print(f"MIC rows                      : {len(mic)}")
print(f"MIC seqs absent from candidates: {len(micseqs-cand_seqs)}")
print(f"candidate peptides lacking MIC : {len(cand_seqs-micseqs)} / {len(cand_seqs)}")
print(f"duplicated MIC index entries   : {mic.index.duplicated().sum()}")
print(f"\ntargets in panel: {mic.shape[1]}")
