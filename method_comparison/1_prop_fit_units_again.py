#%%

# compare methods by looking at the proportion of units that are fit

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import statsmodels.api as sm

plt.style.use('seaborn')

#%% human ac first

# Supplemental Figure 1A

data = pd.read_csv('processed_data_ac.csv')

human_data = data[data.species == 'human']

human_data_filt = human_data[human_data.tau < 1000]
human_data_filt = human_data_filt[human_data.r2 > 0.5]
human_data_filt = human_data_filt[human_data.tau > 10]

# do the plot below of proportion units fit

plt.figure(figsize=(2,2.5))

f = sns.countplot(data=human_data,x='brain_region',hue='species',alpha=0.2,order=['Hippocampus','OFC','Amygdala','ACC'])
g = sns.countplot(data=human_data_filt,x='brain_region',hue='species',order=['Hippocampus','OFC','Amygdala','ACC'])

handles, labels = g.get_legend_handles_labels()
plt.legend(handles[2:], labels[2:],loc='upper right',prop={'size':7})
plt.xlabel('')
plt.ylabel('number of neurons',fontsize=7)
plt.xticks(fontsize=7,rotation=45,ha='right')
plt.yticks(fontsize=7)

plt.title('AC method - Minxha dataset',fontsize=8)

plt.savefig('/Users/zachz/Desktop/humanpropfit.svg',dpi=300)

plt.show()


#%%

titles = ['AC - ITI only','ISI - all spikes','ISI - ITI only','ISI - task only']

# loop over datasets and plot proportion fit neurons
# Supplemental Figure 1B

all_merged_counts = []

for idx, dataset in enumerate(['ac','full','rest','task']):
    
    data_name = 'processed_data_'+dataset+'.csv'

    data = pd.read_csv(data_name)

    data['species'] = pd.Categorical(data['species'], categories = ['mouse','monkey'] , ordered = True)

    # remove human data
    
    data = data[data.species != 'human']
    
    # filter for 10 < tau < 1000 and r2 > 0.5

    filt_data = data[np.logical_and(data.tau < 1000, data.r2 > 0.5)]
    filt_data = filt_data[filt_data.tau > 10]
    
    # below makes a version of Figure 1B for each fit method (without human data)

    # plt.figure(figsize=(2.75,2.75))

    # f = sns.countplot(data=data,x='brain_region',hue='species',alpha=0.2,order=['Hippocampus','OFC','Amygdala','mPFC','ACC'])
    # g = sns.countplot(data=filt_data,x='brain_region',hue='species',order=['Hippocampus','OFC','Amygdala','mPFC','ACC'])

    # handles, labels = g.get_legend_handles_labels()
    # plt.legend(handles[2:], labels[2:],loc='upper right',prop={'size':7})
    # plt.xlabel('')
    # plt.ylabel('number of neurons',fontsize=7)
    # plt.xticks(fontsize=7,rotation=45)
    # plt.yticks(fontsize=7)
    
    # plt.title(titles[idx],fontsize=8)

    # plt.show()
    
    
    # Calculate the total number of neurons for each species and brain area in the original data
    total_counts = data.groupby(['species','brain_region']).size().reset_index(name='total_count')

    # Calculate the number of neurons that were fit for each species and brain area in the filtered data
    fit_counts = filt_data.groupby(['species','brain_region']).size().reset_index(name='fit_count')

    # Merge the total counts and fit counts dataframes
    merged_counts = pd.merge(total_counts, fit_counts, on=['species','brain_region'], how='left')

    # Fill NaN values with 0 (in case there are brain regions with no fit neurons)
    merged_counts['fit_count'] = merged_counts['fit_count'].fillna(0)

    # Calculate the proportion of neurons that were fit
    merged_counts['proportion_fit'] = merged_counts['fit_count'] / merged_counts['total_count']

    # Add method column
    merged_counts['method'] = titles[idx]

    # Append to all_merged_counts list
    all_merged_counts.append(merged_counts)
    
    # Display the result
    # print(merged_counts)
    
    total_total = np.sum(merged_counts.total_count)
    total_fit = np.sum(merged_counts.fit_count)
    
   # print('total proportion fit:',total_fit/total_total)
    
# Combine all merged counts into a single DataFrame
all_merged_counts_df = pd.concat(all_merged_counts)

#%%

# Supplemental Figure 1B

# set plot order

all_merged_counts_df['method'] = pd.Categorical(all_merged_counts_df['method'], categories=['AC - ITI only','ISI - all spikes','ISI - task only','ISI - ITI only'], ordered=True)

# Plot the proportion fit neurons for each species and method

plt.figure(figsize=(3.5,2))
sns.barplot(data=all_merged_counts_df, x='method', y='proportion_fit', hue='species', ci=None)
plt.ylabel('proportion of neurons', fontsize=7)
plt.xlabel('')
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.legend(prop={'size':7})

plt.savefig('/Users/zachz/Desktop/propfit.svg',dpi=300)

plt.show()

#%%