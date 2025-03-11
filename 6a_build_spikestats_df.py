#%% 

# This script calculates various statistics about spike trains, including:
# - number of spikes
# - mean firing rate
# - mean and variance of z-scored firing rate in 1s sliding windows
# - Fano factor
# - burstiness index


# Imports, get filenames
import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import scipy.io as spio

import bottleneck as bn

import mat73

# burstiness - based on Goh and Barabasi, Europhysics Letters 2008

def burstiness_index(spike_times):
    isis = np.diff(spike_times)
    
    # keep values below 1s to exclude breaks between tasks, etc.
    # matches timescale estimation window, also
    #isis = isis[isis < 3]
    
    if len(isis) < 2:
        return np.nan  # Not enough spikes to compute BI
    # keep only ISIs less than 1 sec
    isis = isis[isis < 1]
    mu_isi = np.mean(isis)
    sigma_isi = np.std(isis)
    return (sigma_isi - mu_isi) / (sigma_isi + mu_isi)

#plt.style.use('seaborn')

directory = '/Users/zachz/mounts/backup/Zach/timescales dataset final final backup/datasets_revision'
 
datasets = []

for root, dirs, files in os.walk(directory):
    
    for filename in files:
        
        if filename.startswith('ISI'):
            pass
        
        elif filename.endswith('.mat'):
            
            datasets.append(os.path.join(root, filename))
            
del directory, dirs, files, filename, root

#%% for each dataset, calculate stuff about spiketimes

all_spike_stats = []
     
for ds in range(len(datasets)):
    
    try:
        mat = mat73.loadmat(datasets[ds])
    except:
        mat = spio.loadmat(datasets[ds],squeeze_me=True,simplify_cells=True)
    
    directories = datasets[ds].split('/')

    items = directories[-1].split('_')

    dataset = items[0]
    area = items[1].replace('.mat','')

    spikes = mat['spikes']
    # cell_info = mat['cell_info']
    
    print(dataset, area)
    
    # if type(cell_info) == list:
        
    #     cell_info = cell_info[0]
        
    # elif type(cell_info) == np.ndarray:
        
    #     cell_info = {'dataset': cell_info[0][0],'species': cell_info[0][1],'brain_area':cell_info[0][2]}
    
    for cell in range(len(spikes)):
        
        spiketimes = spikes[cell]
        
        if len(spiketimes) == 0:
            
            pass
        
        else:
        
            # 0 align for ease :)
            
            spiketimes = spiketimes - np.min(spiketimes)
            
            n_spikes = len(spiketimes)
            
            # mean_fr over whole timeseries
            
            mean_fr = len(spiketimes)/np.max(spiketimes)
            
            # mean and var of 1sec sliding windows - z-scored
            
            bins = np.arange(0,np.max(spiketimes),step=0.5) # bin at 500ms
            binned_spikes, edges = np.histogram(spiketimes,bins=bins)
            
            slide = bn.move_sum(binned_spikes, window=2) # sliding window over 1sec (2 windows)
            
            z_slide = slide/mean_fr
            
            mean_z = np.mean(z_slide[1:])
            
            var_z = np.var(z_slide[1:])
            
            # fano factor
            
            fano = var_z / mean_z
            
            # burstiness
            
            burst = burstiness_index(spiketimes)
            
            # assemble
            
            directories = datasets[ds].split('/')
            
            items = directories[-1].split('_')
        
            dataset = items[0]
            area = items[1].replace('.mat','')
            
            if dataset in ['chandravadia','minxha']:
                
                species = 'human'
                
            elif dataset in ['fontanier','meg','stoll','stoll2']:
                species = 'monkey'
                
            elif dataset == 'stein':
                species = 'mouse'

            all_spike_stats.append((dataset,species,area,cell,n_spikes,mean_fr,mean_z,var_z,fano,burst))

# combine into df  
      
data = pd.DataFrame(all_spike_stats,columns=['dataset','species','area','unit','n_spikes','mean_fr','mean_norm_fr','var_norm_fr','fano','burst'])

# rename areas to match other dataframes

amyg = data[(data.area=='amygdala') | (data.area=='bla') | (data.area=='amyg') | (data.area=='AMG')]
hc = data[(data.area=='hippocampus') | (data.area=='hippocampus2') | (data.area=='dg') | (data.area=='ca1') | (data.area=='ca2') | (data.area=='ca3') | (data.area=='hc')]
acc = data[(data.area=='mcc') | (data.area=='dACC') | (data.area=='aca') | (data.area=='ACC')]
mpfc = data[(data.area=='scACC') | (data.area=='ila') | (data.area=='pl')]
ofc = data[(data.area=='OFC') | (data.area=='orb') | (data.area=='a11l') | (data.area=='a11m') | (data.area=='a13l') | (data.area=='a13m')]

acc_g = acc.assign(brain_region = 'ACC')
amyg_g = amyg.assign(brain_region = 'Amygdala')
hc_g = hc.assign(brain_region = 'Hippocampus')
mpfc_g = mpfc.assign(brain_region = 'mPFC')
ofc_g = ofc.assign(brain_region = 'OFC')

grouped_data = pd.concat((acc_g,amyg_g,hc_g,mpfc_g,ofc_g))

brain_regions = ['Hippocampus','Amygdala','OFC','mPFC','ACC']

grouped_data['brain_region'] = pd.Categorical(grouped_data['brain_region'], categories = brain_regions , ordered = True)

#%% load in ISI timescale data

fred_data = pd.read_csv('processed_data.csv')

fred_data['species'] = pd.Categorical(fred_data['species'], categories = ['mouse','monkey','human'] , ordered = True)

filt_data = fred_data[np.logical_and(fred_data.tau < 1000, fred_data.r2 > 0.5)]
filt_data = filt_data[filt_data.tau > 10]

brain_regions = ['Hippocampus','Amygdala','OFC','mPFC','ACC']

filt_data['brain_region'] = pd.Categorical(filt_data['brain_region'], categories = brain_regions , ordered = True)

fred_data = filt_data

#%%

# Merge grouped_data and fred_data on dataset, area, and unitID/unit
merged_data = pd.merge(grouped_data, fred_data, left_on=['dataset', 'area', 'unit'], right_on=['name', 'area', 'unitID'], how='left')

# Select relevant columns and rename them if necessary
merged_data = merged_data[['dataset', 'species_x', 'area', 'unit', 'n_spikes', 'mean_fr', 'mean_norm_fr', 'var_norm_fr', 'fano', 'burst', 'tau', 'lat', 'r2']]
merged_data = merged_data.rename(columns={'species_x': 'species'})


#%%

# loop over units and match

new_df = pd.DataFrame()

for dataset in grouped_data.dataset.unique():
    
    this_dataset = grouped_data[grouped_data.dataset==dataset]
    fred_dataset = fred_data[fred_data.name==dataset]
    
    for brain_area in this_dataset.area.unique():
        
        this_area = this_dataset[this_dataset.area==brain_area]
        fred_area = fred_dataset[fred_dataset.area==brain_area]

        for unit in this_area.unit.unique():
            
            this_unit = this_area[this_area.unit==unit]
            fred_unit = fred_area[fred_area.unitID==unit]
            
            try:
            
                tau = fred_unit['tau'].values[0]
                lat = fred_unit['lat'].values[0]
                r2 = fred_unit['r2'].values[0]
                
                this_unit = this_unit.assign(tau=tau)
                this_unit = this_unit.assign(lat=lat)
                this_unit = this_unit.assign(r2=r2)
                
                new_df = pd.concat((new_df,this_unit),ignore_index=True)
                
            except:
                
                pass

new_df.to_csv('spike_stats.csv')
# %%
