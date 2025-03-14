#%%

# this script generates the brain region hierarchy plots for the single-unit data across species

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn')

# load grouped data, filter, set plot order

data = pd.read_csv('processed_data.csv')

data['species'] = pd.Categorical(data['species'], categories = ['mouse','monkey','human'] , ordered = True)

filt_data = data[np.logical_and(data.tau < 1000, data.r2 > 0.5)]
filt_data = filt_data[filt_data.tau > 10]

brain_regions = ['Hippocampus','OFC','Amygdala','mPFC','ACC']

filt_data['brain_region'] = pd.Categorical(filt_data['brain_region'], categories = brain_regions , ordered = True)

#%%

# Figure 1D, 1E - timescale and latency hierarchies

fig, axs = plt.subplots(1,2,figsize=(4.475,2.75))

sns.lineplot(ax=axs[0],data=filt_data,x='brain_region',y='tau',hue='species',ci=95,markers=True,legend=True,estimator=np.mean)
g = sns.pointplot(ax=axs[0],data=filt_data,x='brain_region',y='tau',hue='species',ci=None,scale=0.4,legend=False)

axs[0].set_xlabel(None)
axs[0].set_xticks(range(len(brain_regions)),brain_regions,ha='right',rotation_mode='anchor')
axs[0].tick_params(axis='x', rotation=45,labelsize=7)
axs[0].tick_params(axis='y',labelsize=7)
axs[0].set_ylabel('timescale (ms)',fontsize=7)

handles, labels = g.get_legend_handles_labels()
axs[0].legend(handles[:3], labels[:3],loc='upper right',prop={'size':7})

sns.lineplot(ax=axs[1],data=filt_data,x='brain_region',y='lat',hue='species',ci=95,markers=True,legend=False,estimator=np.mean)
sns.pointplot(ax=axs[1],data=filt_data,x='brain_region',y='lat',hue='species',ci=None,scale=0.4)


axs[1].set_xlabel(None)
axs[1].set_xticks(range(len(brain_regions)),brain_regions,ha='right',rotation_mode='anchor')
axs[1].tick_params(axis='x', rotation=45,labelsize=7)
axs[1].tick_params(axis='y',labelsize=7)
axs[1].set_ylabel('latency (ms)',fontsize=7)

axs[1].get_legend().remove()

plt.tight_layout()

plt.show()

#%%

# Figure 2B - correlation between timescale and latency

plt.figure(figsize=(2.1,2.1))

sns.regplot(data=filt_data,x='tau',y='lat',ci=95,scatter=False,line_kws={'linestyle':'dashed','color':'black','linewidth':0.75})
sns.kdeplot(data=filt_data,x='tau',y='lat',levels=10,fill=True,color='pink')
plt.xlabel('timescale (ms)',fontsize=7)
plt.ylabel('latency (ms)',fontsize=7)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)

plt.show()
# %%
