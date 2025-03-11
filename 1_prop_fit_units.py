#%%

# this script calculates the proportion of neurons that were successfully fit for each species and brain area

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import statsmodels.api as sm

plt.style.use('seaborn')

# load grouped data

data = pd.read_csv('processed_data.csv')

# order for plotting convenience

data['species'] = pd.Categorical(data['species'], categories = ['mouse','monkey','human'] , ordered = True)

# filter for 10 < tau < 1000 and r2 > 0.5

filt_data = data[np.logical_and(data.tau < 1000, data.r2 > 0.5)]
filt_data = filt_data[filt_data.tau > 10]

#%%

# Figure 1B - proportion neurons fit

plt.figure(figsize=(2.75,2.75))

f = sns.countplot(data=data,x='brain_region',hue='species',alpha=0.2,order=['Hippocampus','OFC','Amygdala','mPFC','ACC']) # total neurons
g = sns.countplot(data=filt_data,x='brain_region',hue='species',order=['Hippocampus','OFC','Amygdala','mPFC','ACC']) # fit neurons

handles, labels = g.get_legend_handles_labels()
plt.legend(handles[3:], labels[3:],loc='upper left',prop={'size':7})
plt.xlabel('')
plt.ylabel('number of neurons',fontsize=7)
plt.xticks(fontsize=7,rotation=45)
plt.yticks(fontsize=7)

plt.show()

#%%

# Calculate the total number of neurons for each species and brain area in the original data
total_counts = data.groupby(['name', 'brain_region']).size().reset_index(name='total_count')

# Calculate the number of neurons that were fit for each species and brain area in the filtered data
fit_counts = filt_data.groupby(['name', 'brain_region']).size().reset_index(name='fit_count')

# Merge the total counts and fit counts dataframes
merged_counts = pd.merge(total_counts, fit_counts, on=['name', 'brain_region'], how='left')

# Fill NaN values with 0 (in case there are brain regions with no fit neurons)
merged_counts['fit_count'] = merged_counts['fit_count'].fillna(0)

# Calculate the proportion of neurons that were fit
merged_counts['proportion_fit'] = merged_counts['fit_count'] / merged_counts['total_count']

# Display the result
print(merged_counts)
