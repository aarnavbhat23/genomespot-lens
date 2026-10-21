"""Build the taxonomy-corrected real ecology extension with executable fits."""
from pathlib import Path
import textwrap
import nbformat as nb
ROOT=Path(__file__).resolve().parents[1]
def md(s):return nb.v4.new_markdown_cell(textwrap.dedent(s).strip())
def code(s):return nb.v4.new_code_cell(textwrap.dedent(s).strip())
cells=[md('''
# Which habitat associations survive stronger controls? — REAL DATA

This exploratory extension corrects an important limitation of notebook 05 and
adds direct habitat pairs, quality restriction, presence and influence checks.
The original supplement has **20 genomes with only `s__` as species**, across six
genera. That unresolved rank cannot establish conspecificity. Notebook 05 pooled
it as a species; those historical fits remain visible but are superseded here.
We exclude unresolved assignments from species-based inference and preserve them
in the 1,964-genome descriptive registry. We do not guess ANI clusters by genus.

The protocol and correction were committed before their respective fits. This
is a secondary exploratory analysis of Macrel small-ORF candidates, a distinct
endpoint from archaeome encrypted APEX windows. Species membership and genomic
GC are incomplete ancestry/composition controls; screened-ORF denominators,
strain phylogeny and source-study controls are still missing. Ecological lineage
sorting can itself be real biology; within-species comparisons answer a narrower
question than total habitat differences.
'''),code('''
from pathlib import Path
import json,sys,tempfile,subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from IPython.display import display
get_ipython().run_line_magic('matplotlib','inline')
ROOT=next(p for p in [Path.cwd(),*Path.cwd().parents] if (p/'src/genomespot_lens').exists())
DATA=ROOT/'examples/ecology-robustness'
BASE=ROOT/'examples/lactobacillaceae-ecology'
provenance=json.loads((DATA/'provenance.json').read_text());assert provenance['data_kind']=='real'
registry=pd.read_csv(DATA/'registry-with-taxonomy-status.csv')
excluded=pd.read_csv(DATA/'unresolved-species-exclusions.csv')
assert len(excluded)==20 and excluded.genus.nunique()==6
assert not excluded.resolved_joint_model_eligible.any()
display(excluded.groupby(['genus','habitat']).size().rename('unresolved genomes'))
print('Descriptive registry:',len(registry),'genomes; resolved labels:',registry[registry.species_resolved].species.nunique())
'''),md('''
## Refit the corrected models and all influence checks

The following cell reruns the real analysis from the frozen registry, including
all 12 direct specifications, four corrected joint models, 35 joint species
omissions and every direct-count species omission. Compare numerical results
with the exported artifacts using floating-point tolerances. No simulated
fallback or remote model service is used.
'''),code('''
with tempfile.TemporaryDirectory() as folder:
    temp=Path(folder)
    for script,out in [('analyze_ecology_robustness.py',temp/'main'),
                       ('audit_ecology_pair_influence.py',temp/'pairs')]:
        completed=subprocess.run([sys.executable,str(ROOT/'scripts'/script),'--data',str(BASE),'--out',str(out)],
            capture_output=True,text=True,check=True)
    effects=pd.read_csv(temp/'main/direct-effects.csv')
    corrected=pd.read_csv(temp/'main/corrected-joint-effects.csv')
    influence=pd.read_csv(temp/'pairs/direct-pair-influence.csv')
    for actual,saved,keys in [
        (effects,DATA/'direct-effects.csv',['effect','ci_low','ci_high','p_value']),
        (corrected,DATA/'corrected-joint-effects.csv',['rate_ratio','ci_low','ci_high']),
        (influence,DATA/'pair-influence/direct-pair-influence.csv',['effect','ci_low','ci_high'])]:
        stored=pd.read_csv(saved)
        assert len(actual)==len(stored)
        for name in keys:assert np.allclose(actual[name],stored[name],rtol=1e-6,atol=1e-9,equal_nan=True)
print('All recomputed estimates agree with frozen real-data artifacts.')
display(corrected[corrected.contrast.eq('fermented food / human')][
    ['model','genomes','species_clusters','rate_ratio','ci_low','ci_high']])
'''),md('''
## Direct pair support and distinct endpoints

For every pair, retain only species represented in both habitats, recalculating
support after completeness >=90% and contamination <5% filtering. Require at least
ten informative species clusters and a full-rank model. Count models retain
individual zeros but exclude all-zero species. Binomial presence models require
within-species variation in presence; their odds ratios use different cohorts
and must not be called count rate ratios. Nonconverged/extreme fits are withheld.

All models include species effects, log genome size, completeness, contamination
and GC. Species-cluster sandwich t intervals use clusters-1 df. Holm correction
covers three contrasts within each endpoint/quality family, conservatively
including unestimable contrasts as p=1. These are exploratory tests.
'''),code('''
display(effects[['endpoint','quality','contrast','genomes','species_clusters','status',
                 'effect','ci_low','ci_high','p_value','p_holm_three_contrasts']])
FIG=ROOT/'study/ecology-comparison/figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'figure.dpi':120,'axes.spines.top':False,'axes.spines.right':False,'font.size':10})
fig,axes=plt.subplots(1,2,figsize=(11,4.6),sharey=True)
contrasts=['human / animal','fermented food / animal','fermented food / human']
for ax,endpoint in zip(axes,['count','presence']):
    for j,contrast in enumerate(contrasts):
        for quality,shift,color in [('published',-.13,'#667684'),('high',.13,'#007f86')]:
            r=effects[(effects.contrast==contrast)&(effects.endpoint==endpoint)&(effects.quality==quality)].iloc[0]
            if r.status=='fit_pass':
                ax.plot([r.ci_low,r.ci_high],[j+shift]*2,color=color,lw=2)
                ax.scatter(r.effect,j+shift,color=color,s=40,label=quality if j==0 else None)
            else:ax.text(1.03,j+shift,'withheld',color=color,va='center',fontsize=8)
    ax.axvline(1,color='black',ls='--',lw=1);ax.set_xscale('log')
    ax.set_xticks([.2,.5,1,2,4],['0.2','0.5','1','2','4']);ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set(xlim=(.12,5),xlabel=('Count rate ratio' if endpoint=='count' else 'Presence odds ratio')+' (95% cluster interval)',
           title='Predicted candidate count' if endpoint=='count' else 'Any predicted candidate')
axes[0].set_yticks(range(3),contrasts);axes[0].invert_yaxis();axes[0].legend(title='Genome quality')
fig.suptitle('REAL DATA: resolved species, direct habitat comparisons',y=1.02)
fig.tight_layout()
for fmt in ['png','svg']:fig.savefig(FIG/f'direct_habitat_robustness.{fmt}',dpi=180,bbox_inches='tight')
plt.show()
'''),md('''
## Influence and the surviving exploratory lead

The high-quality human/animal count contrast is **0.56** (95% cluster interval
0.40–0.77; Holm p≈0.0055 within its three-contrast family), from **379 genomes /
16 resolved species**. Omitting each species produces point estimates 0.53–0.63;
15 of 16 unadjusted omission intervals remain below 1. One crosses 1, so interval
exclusion is not invariant. No omitted species is selected as a new discovery.

Presence is less precise: OR 0.44 (0.15–1.28), with 365 genomes / 13 species.
This does not prove a count-specific biological mechanism; the endpoints have
different cohorts, uncertainty and estimands. The stricter food/animal count
comparison has only nine species, and its presence comparison six, so both are
withheld. The high-quality food/human presence fit is also withheld for extreme
coefficients/nonconvergence. Do not silently report these as null results.
'''),code('''
summary=pd.read_csv(DATA/'pair-influence/pair-influence-summary.csv')
display(summary)
lead=influence[(influence.quality=='high')&(influence.contrast=='human / animal')]
assert len(lead)==16 and lead.status.eq('fit_pass').all()
display(lead[['excluded_species','genomes','effect','ci_low','ci_high']])
print('An omission whose interval reaches 1:')
display(lead[lead.ci_high.ge(1)][['excluded_species','effect','ci_low','ci_high']])
'''),md('''
## Research decision

The corrected joint food/human count estimate is **1.02 (0.75–1.39)**, and the
direct pair estimate **0.92 (0.62–1.36)**. The earlier crude food advantage still
fails these within-species checks. The animal/human count difference is a useful
lead for independent verification, not evidence that a niche increases measured
antimicrobial activity. Full strain phylogeny, proteome composition, exact ORF
screening opportunities and source-study independence remain necessary before
claiming survival of *all* ancestry/composition/coverage controls.

For AllTheBacteria, the deeper public search has not located the original score /
parent / genome export. A released Bakta proteome independently contains 33.3
million distinct 8–50-aa windows, versus 3.9 million reported for the screen.
Its original screen membership is unknown; this is a conditional scope audit,
not proof the paper's total is wrong or a recovered screen. The exact cohort,
fragment-generation procedure and conflicting MIC selection rules require
clarification. See `study/atb-export-search-2026-10-21.md` and the prepared data
request. We do not invent screened zeros or rerun a substitute screen as if it
were the original one.
''')]
book=nb.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'},'language_info':{'name':'python','version':'3.11'}})
nb.write(book,ROOT/'notebooks/06_real_ecology_robustness.ipynb')
print('Built taxonomy-corrected real ecology notebook')
