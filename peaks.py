import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sci
from statannot import add_stat_annotation
import functions as f

sns.set_style("white")
sns.set_context("talk")

# Load the data
dat = pd.read_csv('dataset_Kalifabougou.csv',index_col=0)

# Replace May12 and May13 timepoints with may_timepoint
may_timepoint = 22
dat = dat.replace('May13', may_timepoint, regex=True).replace('May12', may_timepoint, regex=True)

# Convert Timepoint to numeric units of 1 week
t_units = 2
dat['Timepoint'] = pd.to_numeric(dat['Timepoint']) * t_units

# Filter out May12 and May13 timepoints
filtered_dat = dat[~dat['Timepoint'].isin([may_timepoint * t_units])].copy()


#FOLD_CHANGE 

# NA values are filled with the mean of the qPCR values for each kid. qPCR values lower than 1 are set to 1.
first_qpcr = filtered_dat.groupby(['Timepoint', 'Kid']).first().reset_index()
first_qpcr['qPCR'] = first_qpcr['qPCR'].fillna(first_qpcr.groupby('Kid')['qPCR'].transform('mean'))
first_qpcr['qPCR'] = first_qpcr['qPCR'].clip(lower=1)

first_qpcr['FoldChange'] = first_qpcr.groupby('Kid')['qPCR'].pct_change(fill_method=None)+1

# Log transform qPCR counts for more normal distribution
first_qpcr['qPCR'] = np.log10(first_qpcr['qPCR'])


#Calculate the moving avarage of Parasitemia for each kid
first_qpcr['MovingAverage'] = first_qpcr.groupby('Kid')['qPCR'].transform(lambda x: f.rollavg_convolve(x, 3))

# Exclude the first timepoint for each kid (to avoid division by zero)
first_qpcr2= first_qpcr.copy()
first_qpcr= first_qpcr.groupby('Kid').apply(lambda x: x.iloc[1:],include_groups=True ).reset_index(drop=True)
first_qpcr['FoldChange'] = np.log2(first_qpcr['FoldChange'])

#Log2FoldChange plot

plt.figure(figsize=(10, 6))
plt.hist(first_qpcr['FoldChange'], bins=30, edgecolor='black')
plt.xlabel('Fold Change (log2)')
plt.ylabel('Frequency')
plt.show()

#Mean FoldChange per kid plot 
mean_fold_change = first_qpcr.groupby('Kid')['FoldChange'].mean()
plt.figure(figsize=(10, 6))
plt.hist(mean_fold_change, bins=30, edgecolor='black')
plt.axvline(mean_fold_change.mean(), color='red', linestyle='dotted', linewidth=2, label=f'Mean: {mean_fold_change.mean():.2f}')
plt.xlabel('Mean Fold Change (log2)')
plt.ylabel('Frequency')

plt.show()

filtered_dat = filtered_dat.merge(first_qpcr[['Kid', 'Timepoint', 'FoldChange','MovingAverage','qPCR']], on=['Kid', 'Timepoint'], how='left')



# Randomly choose 10 unique kids
random_kids = first_qpcr['Kid'].drop_duplicates().sample(n=10, random_state=42)

# Filter the data for the randomly chosen kids
random_kids_data = first_qpcr[first_qpcr['Kid'].isin(random_kids)]

# Histogram for the randomly chosen kids
plt.figure(figsize=(12, 6))
for kid in random_kids:
    kid_data = random_kids_data[random_kids_data['Kid'] == kid]
    plt.plot(kid_data['Timepoint'], kid_data['qPCR'], marker='o', linestyle='-', label=f'Kid {kid}')

plt.xlabel('Timepoint (weeks)')
plt.ylabel('Parasitemia (Log-transformed qPCR)')
plt.legend(title='Child', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()


# Plot the histogram for the randomly chosen kids
plt.figure(figsize=(12, 6))
for kid in random_kids:
    kid_data = random_kids_data[random_kids_data['Kid'] == kid]
    plt.plot(kid_data['Timepoint'], kid_data['FoldChange'], marker='o', linestyle='-', label=f'Kid {kid}')

plt.xlabel('Timepoint (weeks)')
plt.ylabel('Fold Change')
plt.legend(title='Child', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()

#Average and Standard Deviation of FoldChange at Each Timepoint
timepoints = sorted(first_qpcr['Timepoint'].unique())

fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 12), sharex=True)

palette = sns.color_palette("husl", len(timepoints))
sns.boxplot(data=first_qpcr, x='Timepoint', y='FoldChange', hue='Timepoint', palette=palette, ax=axes[0],legend=False)
axes[0].set_xlabel('')  # Remove x-axis label for the first plot
axes[0].set_ylabel('Log2 Fold Change')

sns.boxplot(data=first_qpcr, x='Timepoint', y='qPCR', hue='Timepoint', palette=palette, ax=axes[1],legend=False)
axes[1].set_xlabel('Timepoint (weeks)')
axes[1].set_ylabel('Parasitemia (Log-transformed qPCR)')

axes[1].set_xticks(range(len(timepoints)))
axes[1].set_xticklabels(timepoints)

plt.tight_layout()
plt.show()




##PEAK DETECTION


# Peak detection derivative-based method
# peak_data_dxqpcr = f.find_peaks_dx(first_qpcr2,'qPCR')
# peak_data_dxfc = f.find_peaks_dx(first_qpcr,'FoldChange')

# Peak detection threshold-based method

peaks_lm_qpcr= f.find_peaks_lm(first_qpcr2,'qPCR',2)
peaks_lm_qpcr2= f.find_peaks_lm(first_qpcr2,'qPCR',2,2)

peaks_lm_fc= f.find_peaks_lm(first_qpcr,'FoldChange',0.02)
peaks_lm_fc2= f.find_peaks_lm(first_qpcr,'FoldChange',0.02,2)

f.plot_peaks_for_random_kids(first_qpcr2, peaks_lm_qpcr2, 'qPCR')

#Peak detection Persistant homology-based method
peak_data_toq1,peak_data_toq2= f.find_peaks_to(first_qpcr2,'qPCR',2)
peak_data_tofc3,peak_data_tofc4= f.find_peaks_to(first_qpcr,'FoldChange',0.02)


