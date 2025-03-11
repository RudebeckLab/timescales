#%%

# this script processes the outputs of the method comparison leading to Figure S1

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import scipy.io as spio
from scipy.io import loadmat

plt.style.use('seaborn')

# load datasets, assign brain regions, set plot order

murray = pd.read_csv('/Users/zachz/Documents/timescales_analysis/2023/paper:poster figures/revision_datasets/ISI_timescales_AC.csv')
isi = pd.read_csv('/Users/zachz/Documents/timescales_analysis/2023/paper:poster figures/revision_datasets/ISI_timescales_full.csv')
isi_rest = pd.read_csv('/Users/zachz/Documents/timescales_analysis/2023/paper:poster figures/revision_datasets/ISI_timescales_rest.csv')
isi_task = pd.read_csv('/Users/zachz/Documents/timescales_analysis/2023/paper:poster figures/revision_datasets/ISI_timescales_task.csv')

def process_data(data):
    amyg = data[(data.area == 'amygdala') | (data.area == 'bla') | (data.area == 'amyg') | (data.area == 'AMG')]
    hc = data[(data.area == 'hippocampus') | (data.area == 'hippocampus2') | (data.area == 'dg') | (data.area == 'ca1') | (data.area == 'ca2') | (data.area == 'ca3') | (data.area == 'hc')]
    acc = data[(data.area == 'mcc') | (data.area == 'dACC') | (data.area == 'aca') | (data.area == 'ACC')]
    mpfc = data[(data.area == 'scACC') | (data.area == 'ila') | (data.area == 'pl')]
    ofc = data[(data.area == 'OFC') | (data.area == 'orb') | (data.area == 'a11l') | (data.area == 'a11m') | (data.area == 'a13l') | (data.area == 'a13m')]

    acc_g = acc.assign(brain_region='ACC')
    amyg_g = amyg.assign(brain_region='Amygdala')
    hc_g = hc.assign(brain_region='Hippocampus')
    mpfc_g = mpfc.assign(brain_region='mPFC')
    ofc_g = ofc.assign(brain_region='OFC')

    grouped_data = pd.concat((acc_g, amyg_g, hc_g, mpfc_g, ofc_g))

    brain_regions = ['Hippocampus', 'Amygdala', 'OFC', 'mPFC', 'ACC']
    grouped_data['brain_region'] = pd.Categorical(grouped_data['brain_region'], categories=brain_regions, ordered=True)
    
    return grouped_data

murray_proc = process_data(murray)
isi_proc = process_data(isi)
isi_rest_proc = process_data(isi_rest)
isi_task_proc = process_data(isi_task)

# save processed outputs

murray_proc.to_csv('/Users/zachz/Documents/timescales_analysis/2023/paper:poster figures/revision_datasets/processed_data_ac.csv', index=False)
isi_proc.to_csv('/Users/zachz/Documents/timescales_analysis/2023/paper:poster figures/revision_datasets/processed_data_full.csv', index=False)
isi_rest_proc.to_csv('/Users/zachz/Documents/timescales_analysis/2023/paper:poster figures/revision_datasets/processed_data_rest.csv', index=False)
isi_task_proc.to_csv('/Users/zachz/Documents/timescales_analysis/2023/paper:poster figures/revision_datasets/processed_data_task.csv', index=False)

#%%

