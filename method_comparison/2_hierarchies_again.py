#%%

# compare hierarchies from outputs of different methods
# not shown in the manuscript or supplement

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

plt.style.use('seaborn')

titles = ['AC - ITI only','ISI - all spikes','ISI - ITI only','ISI - task only']

for idx, dataset in enumerate(['ac','full','rest','task']):
    
    data_name = 'processed_data_'+dataset+'.csv'

    data = pd.read_csv(data_name)

    data['species'] = pd.Categorical(data['species'], categories = ['mouse','monkey','human'] , ordered = True)

    filt_data = data[np.logical_and(data.tau < 1000, data.r2 > 0.5)]
    filt_data = filt_data[filt_data.tau > 10]

    brain_regions = ['Hippocampus','OFC','Amygdala','mPFC','ACC']

    filt_data['brain_region'] = pd.Categorical(filt_data['brain_region'], categories = brain_regions , ordered = True)

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
    
    axs[0].get_legend().remove()
    
    plt.suptitle(titles[idx],fontsize=8)

    plt.tight_layout()

    plt.show()
# %%
