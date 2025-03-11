#%%

# this script plots the distributions of timescales and latencies across brain regions and species

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns

import matplotlib

import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

from scipy.stats import levene

plt.style.use('seaborn')

plt.rcParams['font.size'] = '7'

#%%

# load, filter, plot order

data = pd.read_csv('processed_data.csv')

data['species'] = pd.Categorical(data['species'], categories = ['mouse','monkey','human'] , ordered = True)


filt_data = data[np.logical_and(data.tau < 1000, data.r2 > 0.5)]
filt_data = filt_data[filt_data.tau > 10]

brain_regions = ['Hippocampus','OFC','Amygdala','mPFC','ACC']

filt_data['brain_region'] = pd.Categorical(filt_data['brain_region'], categories = brain_regions , ordered = True)


#%%

# Figure 2C - timescale distributions

fig,axs = plt.subplots(1,5,figsize=(7,2),sharex=True,sharey=True)

for region, ax in zip(brain_regions,axs.ravel()):
        
    this_region = filt_data[filt_data.brain_region==region]
    
    g=sns.kdeplot(ax=ax,data=this_region,x='tau',hue='species',fill=False,common_norm=False,alpha=0.7,log_scale=True)
    
    ax.set_title(region,fontsize=8)
    
    ax.set_xlabel('timescale (ms)',fontsize=7)
    ax.set_ylabel('density',fontsize=7)
    
    ax.tick_params(axis='x',labelsize=7)
    ax.tick_params(axis='y',labelsize=7)  
    # ax.set_xticks([10,100,1000])
    # ax.set_xticklabels(['10','100','1000'])  
    #ax.xaxis.set_label_position('top')
    
    # Calculate and plot the mean for each species
    means = this_region.groupby('species')['tau'].mean()
    colors = sns.color_palette(n_colors=len(means))
    
    for i, ((species, mean), color) in enumerate(zip(means.items(), colors)):
        #ax.axvline(mean, color=color, linestyle='--')
        ax.plot(mean, 1.7 - i * 0.05, 'o', color=color, markersize=3)  # Adjust y-coordinate to avoid overlap
    
    
    plt.setp(g.get_legend().get_texts(), fontsize='7') 
    leg = ax.get_legend()
    leg.set_title('')

    if region!='ACC':
        
        ax.get_legend().remove()
    
plt.tight_layout()
plt.show()

#%%

# levene's test for equal variances

for brain_region in filt_data.brain_region.unique():
    
    if brain_region in ['LAI', 'vlPFC']:
        
        pass
    
    elif brain_region == 'mPFC':
        
        this_region = filt_data[filt_data.brain_region==brain_region]
    
        p = levene(this_region[this_region.species=='mouse'].tau,this_region[this_region.species=='monkey'].tau)
    
        print(brain_region)
        print(p)
    else:
    
        this_region = filt_data[filt_data.brain_region==brain_region]
        
        p = levene(this_region[this_region.species=='mouse'].tau,this_region[this_region.species=='monkey'].tau,this_region[this_region.species=='human'].tau)
        
        print(brain_region)
        print(p)

#%%

# Figure 2D - latency distributions

fig,axs = plt.subplots(1,5,figsize=(7,2),sharex=True,sharey=True)

for region, ax in zip(brain_regions,axs.ravel()):
    
    this_region = filt_data[filt_data.brain_region==region]
    
    g=sns.kdeplot(ax=ax,data=this_region,x='lat',hue='species',fill=False,common_norm=False,alpha=0.7)
    
    ax.set_title(region,fontsize=8)
    
    ax.set_xlabel('latency (ms)',fontsize=7)
    ax.set_ylabel('density',fontsize=7)
    
    #ax.xaxis.set_label_position('top')
    
    ax.tick_params(axis='x',labelsize=7)
    ax.tick_params(axis='y',labelsize=7)
    ax.set_xlim(-50,250)
    
    # Calculate and plot the mean for each species
    means = this_region.groupby('species')['lat'].mean()
    colors = sns.color_palette(n_colors=len(means))
    
    for i, ((species, mean), color) in enumerate(zip(means.items(), colors)):
        #ax.axvline(mean, color=color, linestyle='--')
        ax.plot(mean, 0.049 - i * 0.0015, 'o', color=color, markersize=3)  # Adjust y-coordinate to avoid overlap

    plt.setp(g.get_legend().get_texts(), fontsize='7') 
    leg = ax.get_legend()
    leg.set_title('')
    
    if region!='ACC':
        
        ax.get_legend().remove()
    
plt.tight_layout()
plt.show()

#%%

# levene's test for equal variances

for brain_region in filt_data.brain_region.unique():
    
    if brain_region in ['LAI', 'vlPFC']:
        
        pass
    
    elif brain_region == 'mPFC':
        
        this_region = filt_data[filt_data.brain_region==brain_region]
    
        p = levene(this_region[this_region.species=='mouse'].lat,this_region[this_region.species=='monkey'].lat)
    
        print(brain_region)
        print(p)
    else:
    
        this_region = filt_data[filt_data.brain_region==brain_region]
        
        p = levene(this_region[this_region.species=='mouse'].lat,this_region[this_region.species=='monkey'].lat,this_region[this_region.species=='human'].lat)
        
        print(brain_region)
        print(p)

# %% 

# pairwise levene's tests

from scipy.stats import levene

for brain_region in filt_data.brain_region.unique():
    
    this_region = filt_data[filt_data.brain_region==brain_region]
    
    mouse = this_region[this_region.species=='mouse'].tau.to_numpy()
    monkey = this_region[this_region.species=='monkey'].tau.to_numpy()
    
    if brain_region == 'mPFC':
        
        stat, p = levene(mouse,monkey,center='median')
        
        print(brain_region)
        
        print('F(%i,%i) = %.2f' %(1,len(mouse)+len(monkey),stat))
        
        print('tau p-val = %.3f' %p)
        
    elif brain_region in ['LAI','vlPFC']:
        
        pass
        
    else:
        
        human = this_region[this_region.species=='human'].tau.to_numpy()
        
        stat, p = levene(mouse,monkey,human,center='median')
        
        print(brain_region)
        print('F(%i,%i) = %.2f' %(2,len(mouse)+len(monkey)+len(human),stat))
        print('tau p-val = %.3f' %p)
        
#%%

# pairwise levene's tests

for brain_region in filt_data.brain_region.unique():
    
    this_region = filt_data[filt_data.brain_region==brain_region]
    
    mouse = this_region[this_region.species=='mouse'].lat.to_numpy()
    monkey = this_region[this_region.species=='monkey'].lat.to_numpy()
    
    if brain_region == 'mPFC':
        
        stat, p = levene(mouse,monkey,center='median')
        
        print(brain_region)
        print('F(%i,%i) = %.2f' %(1,len(mouse)+len(monkey),stat))
        print('lat p-val = %.3f' %p)
        
    elif brain_region in ['LAI','vlPFC']:
        
        pass
        
    else:
        
        human = this_region[this_region.species=='human'].lat.to_numpy()
        
        stat, p = levene(mouse,monkey,human,center='median')
        
        print(brain_region)
        print('F(%i,%i) = %.2f' %(2,len(mouse)+len(monkey)+len(human),stat))
        print('lat p-val = %.3f' %p)
        
#%%
