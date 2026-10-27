import sys, re, collections
import pandas as pd

faa, micf, acf = sys.argv[1], sys.argv[2], sys.argv[3]

def headers(p):
    return [l[1:].strip() for l in open(p) if l.startswith('>')]

def records(p):
    h=None; out={}
    for l in open(p):
        l=l.strip()
        if l.startswith('>'): h=l[1:]
        elif l and h: out[h]=out.get(h,'')+l
    return out

cand = records(faa)
acfr = records(acf)
print(f"candidate records: {len(cand)}   distinct seqs: {len(set(cand.values()))}")
print(f"ACF representatives: {len(acfr)} distinct seqs: {len(set(acfr.values()))}")

pat = re.compile(r'^(GC[AF]_\d+\.\d+)_(.+?)_genomic_(.+)_(\d+)$')
rows=[]
bad=0
for h in cand:
    m=pat.match(h)
    if not m: bad+=1; continue
    rows.append({'header':h,'assembly':m.group(1),'asmname':m.group(2),
                 'replicon':m.group(3),'protein_index':int(m.group(4)),'seq':cand[h]})
df=pd.DataFrame(rows)
print(f"parsed: {len(df)}  unparsed: {bad}")
print(f"distinct assemblies: {df.assembly.nunique()}")
print(f"distinct replicons : {df.replicon.nunique()}")
print(f"distinct parent proteins (assembly+replicon+idx): {df.groupby(['assembly','replicon','protein_index']).ngroups}")

pg = df.groupby('assembly').size().sort_values(ascending=False)
print("\n=== candidates per assembly ===")
print(pg.describe().to_string())
print("\ntop 10:"); print(pg.head(10).to_string())
print("\ncount distribution (candidates -> n assemblies):")
print(pg.value_counts().sort_index().head(15).to_string())

mic = pd.read_csv(micf, index_col=0)
print(f"\n=== MIC table: {mic.shape[0]} rows x {mic.shape[1]} targets ===")
micseqs=set(mic.index)
candseqs=set(cand.values()); acfseqs=set(acfr.values())
print(f"MIC seqs in candidate set : {len(micseqs & candseqs)} / {len(micseqs)}")
print(f"MIC seqs in ACF repr set  : {len(micseqs & acfseqs)} / {len(micseqs)}")
print(f"ACF repr seqs with MIC    : {len(acfseqs & micseqs)} / {len(acfseqs)}")
print(f"\nmin-across-targets MIC summary (uM):")
print(mic.min(axis=1).describe().to_string())
print(f"\nrows with any target < 80 uM: {(mic.min(axis=1) < 80).sum()}")
print(f"rows with median target <= 80 uM: {(mic.median(axis=1) <= 80).sum()}")
df.to_csv(sys.argv[4], index=False)
print(f"\nwrote parsed table -> {sys.argv[4]}")
