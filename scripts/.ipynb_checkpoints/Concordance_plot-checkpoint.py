# %%
%load_ext autoreload
%autoreload 2
import pandas as pd
import seaborn as sns
sns.set_theme(style="white")
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
# %%
SNV_overlap = pd.read_csv('../data/SNV_s2i_s2r_i2r.csv')
metadata = pd.read_csv('../data/sampleinfo.txt',sep='\t').iloc[:-1].set_index('SAM id')
# %% Filter out low coverage
coverage = []
for ind, data in SNV_overlap.iterrows():
    coverage.append(metadata.loc[data['SAM id'], f'{data["species"]}_coverage'])

SNV_overlap['coverage'] = coverage
SNV_overlap = SNV_overlap[SNV_overlap['coverage'] < 4]

# %% Consensus
consensus = SNV_overlap[['SAM id', 'species', 'consensus', 'consensus_overlap2', 'consensus_overlap3', 'consensus_overlap4', 'consensus_overlap5']].rename(columns=\
{'consensus': 'total',\
'consensus_overlap2': 's2i',\
'consensus_overlap3': 'i2r',\
'consensus_overlap4': 'uncovered',\
'consensus_overlap5': 's2i+i2r+uncovered'})

# %% Rare
rare = SNV_overlap[['SAM id', 'species', 'rare', 'rare_overlap2', 'rare_overlap3', 'rare_overlap4', 'rare_overlap5']].rename(columns=\
{'rare':'total',\
'rare_overlap2': 's2i',\
'rare_overlap3': 'i2r',\
'rare_overlap4': 'uncovered',\
'rare_overlap5': 's2i+i2r+uncovered'})

# %% Total
total = rare+consensus
total[['SAM id', 'species']] = rare[['SAM id', 'species']]

# %% Total variants (consensus > 0.9 and rare < 0.5)
# Change Color
newcolors = np.array([[6/255,1/255,133/255,1],[128/255,24/255,13/255,1]])
newcmp = ListedColormap(newcolors)
total['concordance'] = total['s2i+i2r+uncovered']/total['total']
data = total
ec = data[data['species']=='E_coli']
kp = data[data['species']=='K_pneumoniae']
concordance = pd.concat([ec['s2i+i2r+uncovered']/ec['total'], kp['s2i+i2r+uncovered']/kp['total']],axis=1).rename(columns={0:'E_coli',1:'K_pneumoniae'})
concordance_long = pd.melt(concordance, value_vars=['E_coli', 'K_pneumoniae'], var_name='Species', value_name='Concordance', ignore_index=True)
fig = plt.figure(figsize=(6,5))
ax = fig.add_subplot(111)
sns.boxplot(ax=ax, x='Species', y='Concordance', data=concordance_long, palette=newcolors,linewidth=3)
sns.stripplot(ax=ax, x = 'Species', y = 'Concordance', data = concordance_long, color='#aaaaa9')

ax.set(ylim=[0.8,1.01])
fig.savefig(f'/home/jovyan/GIS/cre/reports/concordance_total_VCF.pdf')

total.describe()
rare['concordance'] = rare['s2i+i2r+uncovered']/rare['total']
rare.describe()

# %%
