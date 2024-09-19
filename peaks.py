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
dat = pd.read_csv('data\dataset_Kalifabougou.csv',index_col=0)

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








#Average and Standard Deviation of FoldChange at Each Timepoint
timepoints = sorted(first_qpcr['Timepoint'].unique())

fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 12), sharex=True)

palette = sns.color_palette('CMRmap', len(timepoints))
sns.boxplot(
    data=first_qpcr, x='Timepoint', y='FoldChange', hue='Timepoint', palette=palette, ax=axes[0], legend=False,
    boxprops=dict(alpha=0.7)  # Adjust the transparency here
)
axes[0].set_xlabel('')  # Remove x-axis label for the first plot
axes[0].set_ylabel('Log2 Fold Change', fontsize=20)

# Boxplot for qPCR
sns.boxplot(
    data=first_qpcr, x='Timepoint', y='qPCR', hue='Timepoint', palette=palette, ax=axes[1], legend=False,
    boxprops=dict(alpha=0.7)  # Adjust the transparency here
)
axes[1].set_xlabel('Timepoint (weeks)', fontsize=20)
axes[1].set_ylabel('Parasitemia (Log-transformed qPCR)',fontsize=20)


plt.tight_layout()
plt.show()




##PEAK DETECTION


# Peak detection derivative-based method
# peak_data_dxqpcr = f.find_peaks_dx(first_qpcr2,'qPCR')
# peak_data_dxfc = f.find_peaks_dx(first_qpcr,'FoldChange')

# Peak detection threshold-based method

#Think how to incorporate fold change in the threshold method for qPCR

peaks_lm_fc= f.find_peaks_lm(first_qpcr,'FoldChange',0.02) #0.02 is the threshold
peaks_lm_fc2= f.find_peaks_lm(first_qpcr,'FoldChange',0.02,2) #0.02 is the threshold and 2 is the window size of standard deviation
peaks_lm_qpcr= f.find_peaks_lm(first_qpcr2,'qPCR',2) #2 is the threshold
peaks_lm_qpcr2= f.find_peaks_lm(first_qpcr2,'qPCR',2,2) #2 is the threshold and 2 is the window size of standard deviation	

f.plot_peaks_for_random_kids(first_qpcr2, peaks_lm_qpcr2, 'qPCR')

#Peak detection Persistant homology-based method
peak_data_toq1,peak_data_toq2= f.find_peaks_to(first_qpcr2,'qPCR',2)
peak_data_tofc3,peak_data_tofc4= f.find_peaks_to(first_qpcr,'FoldChange',0.02)


#GRAPH AND STATISTICS 

peak_data_toq1['Method'] = 'topology'
peak_data_toq2['Method'] = 'topology'
peaks_lm_qpcr['Method'] = 'local'

# peaks detected by methods

peaksto1= peak_data_toq1[peak_data_toq1['peak'] == True]


f1=f.mergedf(peaksto1, peaks_lm_qpcr)
f2=f.mergedf(peak_data_toq2, peaks_lm_qpcr)


f.plot_heatmap(f1)
f.plot_heatmap(f2)


df= f.create_pivot_df(f2)
filtered_df = df[df['identify_by'] == 'topology'].rename(columns={'topology': 'qPCR'})

f.plot_peaks_for_random_kids(first_qpcr2, filtered_df, 'qPCR')



peaks_lm_qpcr['Data'] = 'qPCR'
peaks_lm_qpcr2['Data'] = 'qPCR'
peaks_lm_fc2['Data'] = 'FoldChange'
peaks_lm_fc['Data'] = 'FoldChange' 

peaks_lm_qpcr['Criteria'] = 'SD+M'
peaks_lm_qpcr2['Criteria'] = '2SD+M'
peaks_lm_fc['Criteria'] = 'SD+M'
peaks_lm_fc2['Criteria'] = '2SD+M'

combined_df = pd.concat([peaks_lm_qpcr, peaks_lm_qpcr2,peaks_lm_fc, peaks_lm_fc2])

sns.set_style("ticks")
g = sns.catplot(
    data=combined_df, kind="swarm",
    x="Criteria", y="Threshold", hue="Criteria", col="Data",
    aspect=.5, palette='CMRmap'
)


for ax in g.axes.flat:
    sns.boxplot(
        x='Criteria', y='Threshold', data=combined_df[combined_df['Data'] == ax.get_title().split(' = ')[1]],
        ax=ax, color='gray', fliersize=0, linewidth=1,
        boxprops=dict(facecolor='none', edgecolor='gray'),
        whiskerprops=dict(color='gray'),
        capprops=dict(color='gray'),
        medianprops=dict(color='gray')
    )
    for label in ax.get_xticklabels():
        label.set_rotation(45)
    ax.set_xlabel('')  

g.set_axis_labels("", "Threshold")
g.set_titles("{col_name}")


plt.tight_layout()
plt.show()


