#%%

# this script assesses how similar the outputs of the 4 timescale estimation methods are

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn')

# load all 4 processed datasets

ac = pd.read_csv('processed_data_ac.csv')
full = pd.read_csv('processed_data_full.csv')
rest = pd.read_csv('processed_data_rest.csv')
task = pd.read_csv('processed_data_task.csv')

# filter each dataset for 10 < tau < 1000 and r2 > 0.5

ac_filt = ac[np.logical_and(ac.tau < 1000, ac.r2 > 0.5)]
ac_filt = ac_filt[ac_filt.tau > 10]

full_filt = full[np.logical_and(full.tau < 1000, full.r2 > 0.5)]
full_filt = full_filt[full_filt.tau > 10]

rest_filt = rest[np.logical_and(rest.tau < 1000, rest.r2 > 0.5)]
rest_filt = rest_filt[rest_filt.tau > 10]

task_filt = task[np.logical_and(task.tau < 1000, task.r2 > 0.5)]
task_filt = task_filt[task_filt.tau > 10]

# Rename columns to ensure unique names
def rename_columns(df, suffix):
    return df.rename(columns={col: f"{col}_{suffix}" for col in df.columns if col not in ['name', 'area', 'unitID']})

ac_filt = rename_columns(ac_filt, 'ac')
full_filt = rename_columns(full_filt, 'full')
rest_filt = rename_columns(rest_filt, 'rest')
task_filt = rename_columns(task_filt, 'task')

# Merge datasets on name, area, and unitID
merged_data = ac_filt.merge(full_filt, on=['name', 'area', 'unitID'])
merged_data = merged_data.merge(rest_filt, on=['name', 'area', 'unitID'])
merged_data = merged_data.merge(task_filt, on=['name', 'area', 'unitID'])

# Calculate correlations
tau_corr = merged_data[['tau_ac', 'tau_full', 'tau_rest', 'tau_task']].corr()
lat_corr = merged_data[['lat_ac', 'lat_full', 'lat_rest', 'lat_task']].corr()

tau_corr = tau_corr.to_numpy()
lat_corr = lat_corr.to_numpy()

#%% test for statistical significance of correlations with fdr correction

from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

tau_ac_full = pearsonr(merged_data['tau_ac'], merged_data['tau_full'])
tau_ac_rest = pearsonr(merged_data['tau_ac'], merged_data['tau_rest'])
tau_ac_task = pearsonr(merged_data['tau_ac'], merged_data['tau_task'])
tau_full_rest = pearsonr(merged_data['tau_full'], merged_data['tau_rest'])
tau_full_task = pearsonr(merged_data['tau_full'], merged_data['tau_task'])
tau_rest_task = pearsonr(merged_data['tau_rest'], merged_data['tau_task'])

p_values = [tau_ac_full[1], tau_ac_rest[1], tau_ac_task[1], tau_full_rest[1], tau_full_task[1], tau_rest_task[1]]
reject, p_values_corrected_tau, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

print(f"Tau AC vs Full: r = %.2f, p_fdr = %.2e" %(tau_ac_full[0], p_values_corrected_tau[0]))
print(f"Tau AC vs Rest: r = %.2f, p_fdr = %.2e" %(tau_ac_rest[0], p_values_corrected_tau[1]))
print(f"Tau AC vs Task: r = %.2f, p_fdr = %.2e" %(tau_ac_task[0], p_values_corrected_tau[2]))
print(f"Tau Full vs Rest: r = %.2f, p_fdr = %.2e" %(tau_full_rest[0], p_values_corrected_tau[3]))
print(f"Tau Full vs Task: r = %.2f, p_fdr = %.2e" %(tau_full_task[0], p_values_corrected_tau[4]))
print(f"Tau Rest vs Task: r = %.2f, p_fdr = %.2e" %(tau_rest_task[0], p_values_corrected_tau[5]))

# repeat for lat

lat_ac_full = pearsonr(merged_data['lat_ac'], merged_data['lat_full'])
lat_ac_rest = pearsonr(merged_data['lat_ac'], merged_data['lat_rest'])
lat_ac_task = pearsonr(merged_data['lat_ac'], merged_data['lat_task'])
lat_full_rest = pearsonr(merged_data['lat_full'], merged_data['lat_rest'])
lat_full_task = pearsonr(merged_data['lat_full'], merged_data['lat_task'])
lat_rest_task = pearsonr(merged_data['lat_rest'], merged_data['lat_task'])


p_values = [lat_ac_full[1], lat_ac_rest[1], lat_ac_task[1], lat_full_rest[1], lat_full_task[1], lat_rest_task[1]]
reject, p_values_corrected_lat, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

print(f"Latency AC vs Full: r = %.2f, p_fdr = %.2e" %(lat_ac_full[0], p_values_corrected_lat[0]))
print(f"Latency AC vs Rest: r = %.2f, p_fdr = %.2e" %(lat_ac_rest[0], p_values_corrected_lat[1]))
print(f"Latency AC vs Task: r = %.2f, p_fdr = %.2e" %(lat_ac_task[0], p_values_corrected_lat[2]))
print(f"Latency Full vs Rest: r = %.2f, p_fdr = %.2e" %(lat_full_rest[0], p_values_corrected_lat[3]))
print(f"Latency Full vs Task: r = %.2f, p_fdr = %.2e" %(lat_full_task[0], p_values_corrected_lat[4]))
print(f"Latency Rest vs Task: r = %.2f, p_fdr = %.2e" %(lat_rest_task[0], p_values_corrected_lat[5]))

correct_pvals_tau = multipletests([0,tau_ac_full[1],tau_ac_rest[1],tau_ac_task[1],tau_ac_full[1],0,tau_full_rest[1],tau_full_task[1],tau_ac_rest[1],tau_full_rest[1],0,tau_rest_task[1],tau_ac_task[1],tau_full_task[1],tau_full_rest[1],0], alpha=0.05, method='fdr_bh')
correct_pvals_tau = correct_pvals_tau[1].reshape((4,4))

correct_pvals_lat = multipletests([0,lat_ac_full[1],lat_ac_rest[1],lat_ac_task[1],lat_ac_full[1],0,lat_full_rest[1],lat_full_task[1],lat_ac_rest[1],lat_full_rest[1],0,lat_rest_task[1],lat_ac_task[1],lat_full_task[1],lat_full_rest[1],0], alpha=0.05, method='fdr_bh')
correct_pvals_lat = correct_pvals_lat[1].reshape((4,4))

# Set diagonal to NaN

for z in range(4):
    tau_corr[z,z] = np.nan
    lat_corr[z,z] = np.nan
    correct_pvals_tau[z,z] = np.nan
    correct_pvals_lat[z,z] = np.nan

# Create custom annotations
tau_annot = np.array([['r = {:.2f}\np = {:.3g}'.format(tau_corr[i, j], correct_pvals_tau[i,j]) for j in range(tau_corr.shape[1])] for i in range(tau_corr.shape[0])])
lat_annot = np.array([['r = {:.2f}\np = {:.3g}'.format(lat_corr[i, j], correct_pvals_lat[i,j]) for j in range(lat_corr.shape[1])] for i in range(lat_corr.shape[0])])

# Supplemental Figure 1C - correlation matrix

fig, axs = plt.subplots(1, 2, figsize=(5.5, 2.9), gridspec_kw={'width_ratios': [1, 1]})  # Keep subplots evenly spaced

sns.heatmap(tau_corr, annot=tau_annot, fmt='', cmap='flare', ax=axs[0], cbar=False, annot_kws={"size": 5})
axs[0].set_xticklabels(['AC\nITI only', 'ISI\nall spikes', 'ISI\nITI only', 'ISI\ntask only'], ha='center', fontsize=7)
axs[0].set_yticklabels(['AC\nITI only', 'ISI\nall spikes', 'ISI\nITI only', 'ISI\ntask only'], ha='center', fontsize=7)
axs[0].set_title('timescale', fontsize=7)

heatmap = sns.heatmap(lat_corr, annot=lat_annot, fmt='', cmap='flare', ax=axs[1], cbar=False, annot_kws={"size": 5})
axs[1].set_xticklabels(['AC\nITI only', 'ISI\nall spikes', 'ISI\nITI only', 'ISI\ntask only'], ha='center', fontsize=7)
axs[1].set_yticklabels(['AC\nITI only', 'ISI\nall spikes', 'ISI\nITI only', 'ISI\ntask only'], ha='center', fontsize=7)
axs[1].set_title('latency', fontsize=7)

cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.75])  # [left, bottom, width, height]
cbar = fig.colorbar(heatmap.collections[0], cax=cbar_ax)
cbar.set_label('correlation (r)', fontsize=7)
cbar.set_ticks([0,0.25, 0.5, 0.75, 1])
cbar.set_ticklabels(['0', '0.25', '0.5', '0.75','1'],fontsize=7)
cbar.mappable.set_clim(0, 1)

# Adjust layout
plt.tight_layout(rect=[0, 0, 0.9, 1])  # Leave space for the colorbar
axs[0].grid(False)
axs[1].grid(False)

# Save figure
plt.savefig('/Users/zachz/Desktop/correlation_across_methods.svg', dpi=300)

plt.show()

#%%

# full area-by-area correlation ISI-full vs AC

merged_data['species_full'] = pd.Categorical(merged_data['species_full'], categories = ['mouse','monkey'] , ordered = True)

plt.figure(figsize=(12, 6))
g = sns.lmplot(data=merged_data, x='tau_full', y='tau_ac', hue='species_full', col='brain_region_full', col_wrap=3)
g.set_titles("{col_name}")
g._legend.set_title("species") 
g.set_axis_labels(x_var = 'ISI_full timescale (ms)', y_var = 'AC timescale (ms)')
#plt.savefig('/Users/zachz/Desktop/tau_ac_full.svg',dpi=300)
plt.show()

#%%

# full area-by-area correlation ISI-full vs ISI-rest

plt.figure(figsize=(12,6))
g = sns.lmplot(data=merged_data, x='tau_full', y='tau_rest', hue='species_full', col='brain_region_full', col_wrap=3)
g.set_titles("{col_name}")
g._legend.set_title("species")
g.set_axis_labels(x_var = 'ISI_full timescale (ms)', y_var = 'ISI_rest timescale (ms)')
#plt.savefig('/Users/zachz/Desktop/tau_rest_full.svg',dpi=300)
plt.show()

#%%

# full area-by-area correlation ISI-full vs ISI-task

plt.figure(figsize=(12,6))
g = sns.lmplot(data=merged_data, x='tau_full', y='tau_task', hue='species_full', col='brain_region_full', col_wrap=3)
g.set_titles("{col_name}")
g._legend.set_title("species") 
g.set_axis_labels(x_var = 'ISI_full timescale (ms)', y_var = 'ISI_task timescale (ms)')
#plt.savefig('/Users/zachz/Desktop/tau_task_full.svg',dpi=300)
plt.show()
#%%

# Supplemental Figure 1D - representative examples of area-by-area correlations

# for each comparison above, plot only the ACC data (matched number of units)

import matplotlib.pyplot as plt
import seaborn as sns

# Filter data
acc_data = merged_data[merged_data['brain_region_full'] == 'ACC']

# Create figure with 3 subplots
fig, axes = plt.subplots(1, 3, figsize=(5.5, 2.5), sharex=True, sharey=True)

# Define plot settings
plot_params = [
    ('tau_ac', 'AC - ITI only\ntimescale (ms)'),
    ('tau_task', 'ISI - task only\ntimescale (ms)'),
    ('tau_rest', 'ISI - ITI only\ntimescale (ms)')
]

# Iterate through subplots
for ax, (y_var, y_label) in zip(axes, plot_params):
    # Plot separate regression lines for each species
    for species in acc_data['species_full'].unique():
        species_data = acc_data[acc_data['species_full'] == species]
        sns.regplot(
            data=species_data, x='tau_full', y=y_var, ax=ax, label=species,
            scatter=True, scatter_kws={'s':5}
        )

    ax.set_xlabel('ISI - all spikes\ntimescale (ms)', fontsize=7)
    ax.set_ylabel(y_label, fontsize=7)
    ax.tick_params(axis='both', labelsize=7)

# Adjust layout and add a shared legend
# handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, labels, title='Species', loc='upper right', fontsize=7)
# for ax in axes:
#     ax.get_legend().remove()  # Remove subplot legends

plt.tight_layout()
plt.savefig('/Users/zachz/Desktop/tau_comparison_ACC.svg', dpi=300)
plt.show()
# %%
