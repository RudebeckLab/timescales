#%%

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import statsmodels.api as sm

plt.style.use('seaborn')

data = pd.read_csv('processed_data.csv')

data['species'] = pd.Categorical(data['species'], categories = ['mouse','monkey','human'] , ordered = True)


filt_data = data[np.logical_and(data.tau < 1000, data.r2 > 0.5)]
filt_data = filt_data[filt_data.tau > 10]

# brain_regions = ['Hippocampus','OFC','Amygdala','mPFC','ACC']

# filt_data['brain_region'] = pd.Categorical(filt_data['brain_region'], categories = brain_regions , ordered = True)

#%%

# brain region data

plt.figure(figsize=(2.75,2.75))

f = sns.countplot(data=data,x='brain_region',hue='species',alpha=0.2,order=['Hippocampus','OFC','Amygdala','mPFC','ACC'])
g = sns.countplot(data=filt_data,x='brain_region',hue='species',order=['Hippocampus','OFC','Amygdala','mPFC','ACC'])

handles, labels = g.get_legend_handles_labels()
plt.legend(handles[3:], labels[3:],loc='upper left',prop={'size':7})
plt.xlabel('')
plt.ylabel('number of neurons',fontsize=7)
plt.xticks(fontsize=7,rotation=45)
plt.yticks(fontsize=7)

plt.show()

#%%

#%%

# Calculate the total number of neurons for each species and brain area in the original data
total_counts = data.groupby(['species', 'brain_region']).size().reset_index(name='total_count')

# Calculate the number of neurons that were fit for each species and brain area in the filtered data
fit_counts = filt_data.groupby(['species', 'brain_region']).size().reset_index(name='fit_count')

# Merge the total counts and fit counts dataframes
merged_counts = pd.merge(total_counts, fit_counts, on=['species', 'brain_region'], how='left')

# Fill NaN values with 0 (in case there are brain regions with no fit neurons)
merged_counts['fit_count'] = merged_counts['fit_count'].fillna(0)

# Calculate the proportion of neurons that were fit
merged_counts['proportion_fit'] = merged_counts['fit_count'] / merged_counts['total_count']

# Display the result
print(merged_counts)

#%%

brain_regions = ['Hippocampus', 'OFC', 'Amygdala', 'mPFC','ACC']

from scipy.stats import sem

# Define the decay function
def decay_function(z, tau, A, B):
    return A * (np.exp(-z / tau) + B)

# Filter the data to include only the specified brain regions
filtered_data = filt_data

# Generate a range of z values
z_values = np.linspace(0, 2000, 2000)  # Adjust the range and number of points as needed

# Define the order of brain regions

# get mean latencies
mean_lats = []
        
means = []
sems = []
all_species = []
all_areas  = []

for species in ['mouse','monkey','human']:
    
    this_species = filtered_data[filtered_data['species'] == species]
    
    for region in brain_regions:
        
        if species == 'human' and region == 'mPFC':
            pass
        
        else:
        
            all_decay_curves = []
            all_lats = []
                            
            this_region = this_species[this_species['brain_region'] == region]
            
            for neuron in range(len(this_region)):
                
                tau = this_region['tau'].iloc[neuron]
                A = this_region['A'].iloc[neuron]
                B = this_region['B'].iloc[neuron]
                
                lat = this_region['lat'].iloc[neuron]
                
                decay_curve = decay_function(z_values, tau, A, B)
                
                # nan all values before latency
                decay_curve[:int(lat)] = np.nan #
                
                # delete nans
                decay_curve = decay_curve[~np.isnan(decay_curve)]
                
                # cap at len 800
                decay_curve = decay_curve[:600]
                
                all_decay_curves.append(decay_curve)
                
                all_lats.append(lat)
        
            mean_decay_curve = np.mean(all_decay_curves, axis=0) 
            sem_decay_curve = sem(all_decay_curves, axis=0)
            mean_lats.append(int(np.mean(all_lats)))
            means.append(mean_decay_curve)
            sems.append(sem_decay_curve)
            all_species.append(species)
            all_areas.append(region)

#%%

# plot average decay curves


from scipy.stats import sem

# Define the order of brain regions
brain_regions = ['Hippocampus', 'OFC', 'Amygdala', 'mPFC','ACC']

z_values = z_values[:600] # cap at 800

# Create subplots
fig, axes = plt.subplots(1, len(brain_regions), sharex=True, figsize=(7,2.4))

# Ensure axes is always a list
if len(brain_regions) == 1:
    axes = [axes]

# Iterate over each brain region
for ax, brain_region in zip(axes, brain_regions):
    # Filter the means, species, and areas for the current brain region
    for mean_decay_curve, sem_decay_curve, species, area,mean_lat in zip(means, sems, all_species, all_areas,mean_lats):
        if area == brain_region:
            # Calculate the standard error of the mean            
            # Plot the averaged decay curve
            ax.plot(z_values+int(mean_lat), mean_decay_curve, label=f'{species} ($\\tau$ =  ms)')
            ax.fill_between(z_values+int(mean_lat), mean_decay_curve - sem_decay_curve, mean_decay_curve + sem_decay_curve, alpha=0.3)
    
    ax.set_title(f'{brain_region}', fontsize=7)
    ax.set_xlabel('time lag (ms)', fontsize=7)
    ax.set_ylim(0, 1)
    
    if brain_region == 'Hippocampus':
        ax.set_ylabel('autocorrelation (a.u.)', fontsize=7)
    else:
        ax.set_ylabel('')   
        ax.set_yticklabels([])
    ax.tick_params(axis='both', labelsize=7)
    #ax.legend(fontsize=7)

plt.tight_layout()
plt.show()
                
#%%
