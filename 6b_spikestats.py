#%%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns

import matplotlib

import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
import statsmodels.api as sm


from scipy.stats import zscore

plt.style.use('seaborn-v0_8')
plt.rcParams['font.size'] = '7'

data = pd.read_csv('spike_stats.csv')

data['z_n_spikes'] = data.groupby(['species','brain_region']).n_spikes.transform(zscore, ddof=1)
data['z_mean_fr'] = data.groupby(['species','brain_region']).mean_fr.transform(zscore, ddof=1)
data['z_mean_norm_fr'] = data.groupby(['species','brain_region']).mean_norm_fr.transform(zscore, ddof=1)
data['z_fano'] = data.groupby(['species','brain_region']).fano.transform(zscore, ddof=1)
data['z_prop_burst'] = data.groupby(['species','brain_region']).prop_burst.transform(zscore, ddof=1)
data['z_prop_pause'] = data.groupby(['species','brain_region']).prop_pause.transform(zscore, ddof=1)
#%% get PEV and p-values for each predictor

pvals = []
pevs = []

model = smf.ols('tau ~ species + brain_region',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(np.sum(anova['sum_sq'][:-1])/anova['sum_sq'].sum())

#

model = smf.ols('tau ~ species + brain_region + z_n_spikes',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(anova['sum_sq'][-2]/anova['sum_sq'].sum())

#

model = smf.ols('tau ~ species + brain_region + z_mean_fr',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(anova['sum_sq'][-2]/anova['sum_sq'].sum())

#

model = smf.ols('tau ~ species + brain_region + z_mean_norm_fr',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(anova['sum_sq'][-2]/np.sum(anova['sum_sq']))

#

model = smf.ols('tau ~ species + brain_region + z_fano',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(anova['sum_sq'][-2]/anova['sum_sq'].sum())

#

model = smf.ols('tau ~ species + brain_region + z_prop_burst',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(anova['sum_sq'][-2]/anova['sum_sq'].sum())

#

model = smf.ols('tau ~ species + brain_region + z_prop_pause',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(anova['sum_sq'][-2]/anova['sum_sq'].sum())

pevs = np.array(pevs)*100


#%%

plt.figure(figsize=(4,2))

plt.bar(x=range(len(pevs)),height=pevs)
plt.yscale('log')

plt.title('Timescale',size=7)

plt.ylabel('percent explained variance',fontsize=7)
plt.yticks([0.0001, 0.01, 1,10,100], [0.0001, 0.01, 1,10,100], fontsize=7)

plt.xticks(range(7),['area +\nspecies', 'n spikes','mean\nFR','norm\nFR','fano','prop\nburst','prop\npause'],fontsize=7,rotation=0,ha='center')

# add p-values above each bar

for i in range(len(pvals)):
    plt.text(i,pevs[i],f'p={pvals[i]:.4f}',ha='center',va='bottom',fontsize=5)

plt.show()

#%% repeat above for latency

pvals = []
pevs = []

model = smf.ols('lat ~ species + brain_region',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(np.sum(anova['sum_sq'][:-1])/anova['sum_sq'].sum())

#

model = smf.ols('lat ~ species + brain_region + z_n_spikes',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(anova['sum_sq'][-2]/anova['sum_sq'].sum())

#

model = smf.ols('lat ~ species + brain_region + z_mean_fr',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(anova['sum_sq'][-2]/anova['sum_sq'].sum())

#

model = smf.ols('lat ~ species + brain_region + z_mean_norm_fr',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(anova['sum_sq'][-2]/np.sum(anova['sum_sq']))

#

model = smf.ols('lat ~ species + brain_region + z_fano',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(anova['sum_sq'][-2]/anova['sum_sq'].sum())

#

model = smf.ols('lat ~ species + brain_region + z_prop_burst',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(anova['sum_sq'][-2]/anova['sum_sq'].sum())

#

model = smf.ols('lat ~ species + brain_region + z_prop_pause',data=data)
res = model.fit()

pvals.append(res.pvalues[-1])

anova = anova_lm(res)
anova

pevs.append(anova['sum_sq'][-2]/anova['sum_sq'].sum())

pevs = np.array(pevs)*100

#%%

plt.figure(figsize=(4,2))

plt.bar(x=range(len(pevs)),height=pevs)
plt.yscale('log')

plt.title('Latency',size=7)

plt.ylabel('percent explained variance',fontsize=7)
plt.yticks([0.0001, 0.01, 1,10,100], [0.0001, 0.01, 1,10,100], fontsize=7)

plt.xticks(range(7),['area +\nspecies', 'n spikes','mean\nFR','norm\nFR','fano','prop\nburst','prop\npause'],fontsize=7,rotation=0,ha='center')

# add p-values above each bar

# make the p in p-value italics
for i in range(len(pvals)):
    plt.text(i,pevs[i],f'p={pvals[i]:.4f}',ha='center',va='bottom',fontsize=5)

plt.show()

# %%
