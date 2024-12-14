#Author: Ingrid Vanessa Mora Sanchez 
#Email: i.moras@uniandes.edu.co 
#Portugal Lab 


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sci
from statannot import add_stat_annotation
import seaborn.objects as so
import functions as f
import os
import scipy.stats as stats
# Load the data
data = pd.read_csv('..\\data\\dataset_Kalifabougou.csv')	

# Replace May12 and May13 timepoints with may_timepoint
may_timepoint = 13
data = data.replace('May13', may_timepoint, regex=True).replace('May12', 0, regex=True)
data['Timepoint'] = pd.to_numeric(data['Timepoint']) 

# Filter out May12 and May13 timepoints
filtered_data = data[~data['Timepoint'].isin([may_timepoint,0])].copy()
first_data= filtered_data.drop_duplicates(subset=['Timepoint', 'Kid'], keep='first').copy()


# NA values are filled with the mean of the qPCR values for each kid. qPCR values lower than 1 are set to 1.

first_data['qPCR'] = first_data['qPCR'].fillna(first_data.groupby('Kid')['qPCR'].transform('mean'))
first_data['qPCR'] = first_data['qPCR'].clip(lower=1)

# Log transform qPCR counts for more normal distribution
first_data['log2_qPCR'] = np.log2(first_data['qPCR'])


first_data['log10_qPCR'] = np.log10(first_data['qPCR'])
first_data['log2FoldChange'] = first_data.groupby('Kid')['log2_qPCR'].diff()
first_data['log10FoldChange'] = first_data.groupby('Kid')['log10_qPCR'].diff()

#Calculate the moving avarage of Parasitemia for each kid
first_data['MovingAverage'] = first_data.groupby('Kid')['qPCR'].transform(lambda x: f.rollavg_convolve(x, 3))

# Exclude the first timepoint for each kid (to avoid division by zero)
first_data_ft= first_data.copy() #Dataframe with the first timepoint for each kid
first_data= first_data.groupby('Kid').apply(lambda x: x.iloc[1:],include_groups=True ).reset_index(drop=True)
first_data['FoldChange'] = first_data['log2FoldChange'].apply(lambda x: pow(2, x))



#Mean FoldChange per kid plot 
mean_fold_change = first_data.groupby('Kid')['log2_qPCR'].std()
plt.figure(figsize=(10, 6))
plt.hist(mean_fold_change, bins=30, edgecolor='black')
plt.axvline(mean_fold_change.mean(), color='red', linestyle='dotted', linewidth=2, label=f'Mean: {mean_fold_change.mean():.2f}')
quantiles = mean_fold_change.quantile([0.25, 0.5, 0.75])
for q in quantiles.index:
    plt.axvline(quantiles[q], color='blue', linestyle='dotted', linewidth=2, label=f'{int(q*100)}th percentile: {quantiles[q]:.2f}')
plt.xlabel('Mean Fold Change (log2)')

# Check if the distribution is normal
k2, p = stats.normaltest(mean_fold_change)
if p < 0.05:
    plt.title('Mean Fold Change Distribution (Not Normal)')
else:
    plt.title('Mean Fold Change Distribution (Normal)')
plt.ylabel('Frequency')

plt.show()

filtered_data = filtered_data.merge(first_data[['Kid', 'Timepoint', 'log2FoldChange','MovingAverage','log2_qPCR']], on=['Kid', 'Timepoint'], how='left')


#Average and Standard Deviation of FoldChange at Each Timepoint
timepoints = sorted(first_data['Timepoint'].unique())

fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 12), sharex=True)

palette = sns.color_palette('CMRmap', len(timepoints))
sns.boxplot(
    data=first_data, x='Timepoint', y='log2FoldChange', hue='Timepoint', palette=palette, ax=axes[0], legend=False,
    boxprops=dict(alpha=0.7)  
)
axes[0].set_xlabel('')
axes[0].set_ylabel('Log2 Fold Change', fontsize=20)

# Boxplot for qPCR
sns.boxplot(
    data=first_data, x='Timepoint', y='log2_qPCR', hue='Timepoint', palette=palette, ax=axes[1], legend=False,
    boxprops=dict(alpha=0.7)  
)
axes[1].set_xlabel('Timepoint (weeks)', fontsize=20)
axes[1].set_ylabel('Parasitemia (Log2-transformed qPCR)',fontsize=20)


plt.tight_layout()
plt.show()




##PEAK DETECTION
peaks_lm_qpcr= f.find_peaks_lm(first_data_ft,'log2_qPCR',np.log2(100),1) 
peaks_lm_qpcr2= f.find_peaks_lm(first_data_ft,'log2_qPCR',np.log2(100),2) #2 is the threshold
peaks_lm_qpcrs= peaks_lm_qpcr[peaks_lm_qpcr['peak'] == True]
peaks_lm_qpcr3= f.find_peaks_lm(first_data_ft,'log2_qPCR',np.log2(100),3)
peaks_lm_qpcrs2= peaks_lm_qpcr2[peaks_lm_qpcr2['peak'] == True]
peaks_lm_qpcr4= f.find_peaks_lm(first_data_ft,'log2_qPCR',np.log2(100),4)




#Peak detection Persistant homology-based method
peak_data_toq1,peak_data_toq2= f.find_peaks_to(first_data,'log2_qPCR', np.log2(100))



#GRAPH AND STATISTICS 

peak_data_toq1['Method'] = 'topology'
peak_data_toq2['Method'] = 'topology'
peaks_lm_qpcr['Method'] = 'local'

# peaks detected by methods

peaksto1= peak_data_toq1[peak_data_toq1['peak'] == True]


f1=f.mergedf(peaksto1, peaks_lm_qpcr)
f2=f.mergedf(peak_data_toq2, peaks_lm_qpcr)

mean_f2_kids = pd.crosstab(f2['Kid'], f2['Method']).mean()
mean_f2_tmp = pd.crosstab(f2['Timepoint'], f2['Method']).mean()


f.plot_heatmap(f1)
f.plot_heatmap(f2)


df= f.create_pivot_df(f2)


peaks_lm_qpcr['Data'] = 'log2qPCR'
peaks_lm_qpcr2['Data'] = 'log2qPCR'
peaks_lm_qpcr3['Data'] = 'log2qPCR'
peaks_lm_qpcr4['Data'] = 'log2qPCR'

peaks_lm_qpcr['Criteria'] = 'LOD'
peaks_lm_qpcr2['Criteria'] = 'MEAN'
peaks_lm_qpcr3['Criteria'] = 'MEDIAN'
peaks_lm_qpcr4['Criteria'] = '75P'

combined_df = pd.concat([peaks_lm_qpcr, peaks_lm_qpcr2, peaks_lm_qpcr3,peaks_lm_qpcr4])
count_df = combined_df.groupby('Criteria').size().reset_index(name='count')
sns.set_style("ticks")

palette=[ '#75A589', '#ABE891','#365F6A', '#408F8A']
count_df = combined_df.groupby('Criteria').size().reset_index(name='count')

sns.set_style("ticks")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))

sns.swarmplot(
    data=combined_df, x="Criteria", y='height', hue="Criteria", palette=palette, ax=ax1
)

for i, criteria in enumerate(combined_df['Criteria'].unique()):
    sns.boxplot(
        x='Criteria', y='height', data=combined_df[combined_df['Criteria'] == criteria], ax=ax1,
        color=palette[i], fliersize=0, linewidth=1,
        boxprops=dict(facecolor='none', edgecolor=palette[i]),
        whiskerprops=dict(color=palette[i]),
        capprops=dict(color=palette[i]),
        medianprops=dict(color=palette[i])
    )

sns.barplot(
    data=count_df, x='Criteria', y='count', palette=palette, ax=ax2
)

# Ajustar las propiedades del gráfico de enjambre y caja
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45)
ax1.set_xlabel('')
ax1.set_ylabel('height')

# Ajustar las propiedades del gráfico de barras
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45)
ax2.set_xlabel('Criteria')
ax2.set_ylabel('Number of Peaks')

plt.tight_layout()
plt.show()





kids_peaks = first_data_ft[first_data_ft['Kid'].isin(df['Kid'])]
kids = pd.merge(kids_peaks, df, on=['Kid', 'Timepoint'], how='left')
kids=kids.rename(columns={'log2_qPCR_x':'log2_qPCR'})
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from seaborn import axes_style
so.Plot.config.theme.update(axes_style("ticks"))
so.Plot.config.display["format"] = "svg"

with PdfPages('..\\outcome\\kids_peaks_plots.pdf') as pdf:
    for kid in kids['Kid'].unique():
        kid_data = kids[kids['Kid'] == kid]

        plt.figure(figsize=(10, 5))  
        plot = (
        so.Plot(kid_data, x="Timepoint", y="log2_qPCR")
        .add(so.Line(color=".2"))
        .add(so.Dot(), color="Method", fill='peak',marker='falsetype')
        .label(x="Timepoint (weeks)", y="Log2(qPCR)")
        .limit(x=(2, 24),y=(0,18))  
        )
        fig, ax = plt.subplots(figsize =(2.8,2))
        plot.on(ax).plot()
        ax.set_xticks(range(2, 24, 2))
        
        
        plt.title(f'Kid {kid}')
        
        plt.tight_layout()
        pdf.savefig(dpi = 300, orientation = 'portrait', bbox_inches = 'tight')
        plt.close()

print("kids_peaks_plots.pdf")


sns.FacetGrid(kids, col="Kid", row="Method", hue="peak", margin_titles=True).map(plt.plot, "Timepoint", "log2_qPCR", marker="o").add_legend()






#merged peaks data with filtered_dat

merged = pd.merge(data, peaks_lm_qpcr, on=['Kid', 'Timepoint'], how='left')

merged['peak'] = merged['peak'].fillna(False).astype(int)




merged.to_csv('..\\outcome\\peaks.csv', index=False)


#Stats

#Arreglar esta grafica


merged_1= merged.groupby(['Timepoint', 'Kid']).first().reset_index()
merged_2 = merged_1[(merged_1['peak'] == 1)]

promedios_picos = merged_2.groupby(['Kid', 'Type'])['peak'].count().reset_index()

plt.figure()
sns.histplot(data=promedios_picos, x='peak', hue='Type', stat="frequency", common_norm=False)
plt.xlabel('Average number of peaks per kid')
plt.ylabel('Frequency')
plt.show()