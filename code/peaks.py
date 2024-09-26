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
dat = pd.read_csv('..\\data\\dataset_Kalifabougou.csv',index_col=0)

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

# Log transform qPCR counts for more normal distribution
first_qpcr['log2_qPCR'] = np.log2(first_qpcr['qPCR'])


first_qpcr['log10_qPCR'] = np.log10(first_qpcr['qPCR'])
first_qpcr['log2FoldChange'] = first_qpcr.groupby('Kid')['log2_qPCR'].diff()
first_qpcr['log10FoldChange'] = first_qpcr.groupby('Kid')['log10_qPCR'].diff()
#Calculate the moving avarage of Parasitemia for each kid
first_qpcr['MovingAverage'] = first_qpcr.groupby('Kid')['qPCR'].transform(lambda x: f.rollavg_convolve(x, 3))

# Exclude the first timepoint for each kid (to avoid division by zero)
first_qpcr2= first_qpcr.copy() #Dataframe with the first timepoint for each kid
first_qpcr= first_qpcr.groupby('Kid').apply(lambda x: x.iloc[1:],include_groups=True ).reset_index(drop=True)
first_qpcr['FoldChange'] = first_qpcr['log2FoldChange'].apply(lambda x: pow(2, x))

#Log2FoldChange plot

plt.figure(figsize=(10, 6))
plt.hist(first_qpcr['log2FoldChange'], bins=30, edgecolor='black')
plt.xlabel('Fold Change (log2)')
plt.ylabel('Frequency')
plt.show()

#Mean FoldChange per kid plot 
mean_fold_change = first_qpcr.groupby('Kid')['log2FoldChange'].mean()
plt.figure(figsize=(10, 6))
plt.hist(mean_fold_change, bins=30, edgecolor='black')
plt.axvline(mean_fold_change.mean(), color='red', linestyle='dotted', linewidth=2, label=f'Mean: {mean_fold_change.mean():.2f}')
plt.xlabel('Mean Fold Change (log2)')
plt.ylabel('Frequency')

plt.show()

filtered_dat = filtered_dat.merge(first_qpcr[['Kid', 'Timepoint', 'log2FoldChange','MovingAverage','log2_qPCR']], on=['Kid', 'Timepoint'], how='left')





#Average and Standard Deviation of FoldChange at Each Timepoint
timepoints = sorted(first_qpcr['Timepoint'].unique())

fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 12), sharex=True)

palette = sns.color_palette('CMRmap', len(timepoints))
sns.boxplot(
    data=first_qpcr, x='Timepoint', y='log2FoldChange', hue='Timepoint', palette=palette, ax=axes[0], legend=False,
    boxprops=dict(alpha=0.7)  # Adjust the transparency here
)
axes[0].set_xlabel('')  # Remove x-axis label for the first plot
axes[0].set_ylabel('Log2 Fold Change', fontsize=20)

# Boxplot for qPCR
sns.boxplot(
    data=first_qpcr, x='Timepoint', y='log2_qPCR', hue='Timepoint', palette=palette, ax=axes[1], legend=False,
    boxprops=dict(alpha=0.7)  # Adjust the transparency here
)
axes[1].set_xlabel('Timepoint (weeks)', fontsize=20)
axes[1].set_ylabel('Parasitemia (Log2-transformed qPCR)',fontsize=20)


plt.tight_layout()
plt.show()




##PEAK DETECTION


# Peak detection derivative-based method
# peak_data_dxqpcr = f.find_peaks_dx(first_qpcr2,'qPCR')
# peak_data_dxfc = f.find_peaks_dx(first_qpcr,'FoldChange')

# Peak detection threshold-based method

#Think how to incorporate fold change in the threshold method for qPCR

peaks_lm_qpcr= f.find_peaks_lm(first_qpcr2,'log2_qPCR',np.log2(100),1) 
peaks_lm_qpcr2= f.find_peaks_lm(first_qpcr2,'log2_qPCR',np.log2(100),2) #2 is the threshold
peaks_lm_qpcrs= peaks_lm_qpcr[peaks_lm_qpcr['peak'] == True]
peaks_lm_qpcr3= f.find_peaks_lm(first_qpcr2,'log2_qPCR',np.log2(100),3)
peaks_lm_qpcrs2= peaks_lm_qpcr2[peaks_lm_qpcr2['peak'] == True]
peaks_lm_qpcr4= f.find_peaks_lm(first_qpcr2,'log2_qPCR',np.log2(100),4)


f.plot_peaks_for_random_kids(first_qpcr2, peaks_lm_qpcr, 'log2_qPCR')

#Peak detection Persistant homology-based method
peak_data_toq1,peak_data_toq2= f.find_peaks_to(first_qpcr2,'log2_qPCR', np.log2(100))



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


df= f.create_pivot_df(f1)


peaks_lm_qpcr['Data'] = 'log2qPCR'
peaks_lm_qpcr2['Data'] = 'log2qPCR'
peaks_lm_qpcr3['Data'] = 'log2qPCR'
peaks_lm_qpcr4['Data'] = 'log2qPCR'

peaks_lm_qpcr['Criteria'] = 'LOD'
peaks_lm_qpcr2['Criteria'] = 'MEAN'
peaks_lm_qpcr3['Criteria'] = 'MEDIAN'
peaks_lm_qpcr4['Criteria'] = '75P'

combined_df = pd.concat([peaks_lm_qpcr, peaks_lm_qpcr2, peaks_lm_qpcr3,peaks_lm_qpcr4])

sns.set_style("ticks")
g = sns.catplot(
    data=combined_df, kind="swarm",
    x="Criteria", y='height', hue="Criteria", col="Data",
    aspect=.5, palette='CMRmap'
)


for ax in g.axes.flat:
    sns.boxplot(
        x='Criteria', y='height', data=combined_df[combined_df['Data'] == ax.get_title().split(' = ')[1]],
        ax=ax, color='gray', fliersize=0, linewidth=1,
        boxprops=dict(facecolor='none', edgecolor='gray'),
        whiskerprops=dict(color='gray'),
        capprops=dict(color='gray'),
        medianprops=dict(color='gray')
    )
    for label in ax.get_xticklabels():
        label.set_rotation(45)
    ax.set_xlabel('')  

g.set_axis_labels("", "height")
g.set_titles("{col_name}")


plt.tight_layout()
plt.show()


kids_peaks = first_qpcr2[first_qpcr2['Kid'].isin(df['Kid'])]
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

with PdfPages('..\\outcome\\kids_peaks_plots.pdf') as pdf:
    for kid in df['Kid'].unique():
        kid_data = kids_peaks[kids_peaks['Kid'] == kid]
        kid_df = df[df['Kid'] == kid]

        plt.figure()  
        
        # Scatter plot for 'both'
        both_data = kid_df[kid_df['identify_by'] == 'both']
        if not both_data.empty:
            plt.scatter(both_data[both_data['peak'] == True]['Timepoint'], both_data[both_data['peak'] == True]['log2_qPCR'], color='red', marker='o', label='Both (Peak)', zorder=3)
            plt.scatter(both_data[both_data['peak'] == False]['Timepoint'], both_data[both_data['peak'] == False]['log2_qPCR'], color='red', marker='^', label='Both (Non-Peak)', zorder=3)

        # Scatter plot for 'topology'
        topology_data = kid_df[kid_df['identify_by'] == 'topology']
        if not topology_data.empty:
            plt.scatter(topology_data[topology_data['peak'] == True]['Timepoint'], topology_data[topology_data['peak'] == True]['log2_qPCR'], color='blue', marker='o', label='Topology (Peak)', zorder=3)
            plt.scatter(topology_data[topology_data['peak'] == False]['Timepoint'], topology_data[topology_data['peak'] == False]['log2_qPCR'], color='blue', marker='^', label='Topology (Non-Peak)', zorder=3)

        # Scatter plot for 'local'
        local_data = kid_df[kid_df['identify_by'] == 'local']
        if not local_data.empty:
            plt.scatter(local_data[local_data['peak'] == True]['Timepoint'], local_data[local_data['peak'] == True]['log2_qPCR'], color='green', marker='o', label='Local (Peak)', zorder=3)
            plt.scatter(local_data[local_data['peak'] == False]['Timepoint'], local_data[local_data['peak'] == False]['log2_qPCR'], color='green', marker='^', label='Local (Non-Peak)', zorder=3)

        # Plot the line plot for the kid
        plt.plot(kid_data['Timepoint'], kid_data['log2_qPCR'], marker='o', linestyle='-', color='black', label=f'Kid {kid}', zorder=1)
        
        plt.xlabel('Timepoint (weeks)')
        plt.ylabel('Parasitemia (Log-transformed qPCR)')
        
        # Create custom legend handles
        legend_elements = [
            Line2D([0], [0], color='red', marker='o', linestyle='None', label='Both (Peak)'),
            Line2D([0], [0], color='red', marker='^', linestyle='None', label='Both (Non-Peak)'),
            Line2D([0], [0], color='blue', marker='o', linestyle='None', label='Topology (Peak)'),
            Line2D([0], [0], color='blue', marker='^', linestyle='None', label='Topology (Non-Peak)'),
            Line2D([0], [0], color='green', marker='o', linestyle='None', label='Local (Peak)'),
            Line2D([0], [0], color='green', marker='^', linestyle='None', label='Local (Non-Peak)'),
            Line2D([0], [0], color='black', marker='o', linestyle='-', label=f'Kid {kid}')
        ]
        
        plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.title(f'Kid {kid}')
        
        plt.tight_layout()
        pdf.savefig()
        plt.close()

print("kids_peaks_plots.pdf")


#merged peaks data with filtered_dat
