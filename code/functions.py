
from matplotlib import pyplot as plt
import numpy as np
import scipy.signal as sig
from scipy.signal import find_peaks, peak_prominences
import pandas as pd
import seaborn as sns
import seaborn.objects as so
from findpeaks import findpeaks
from matplotlib.colors import ListedColormap
from matplotlib.backends.backend_pdf import PdfPages

def rollavg_convolve(a, n):
    'scipy.convolve'
    assert n % 2 == 1
    if len(a) < n:
        return np.full_like(a, np.nan) 
    convolved = sig.convolve(a, np.ones(n, dtype='float') / n, 'same')
    pad_size = n // 2
    return np.pad(convolved[pad_size:-pad_size], (pad_size, pad_size), mode='edge')


def find_peaks_lm(df, col, lod, opt=1):
    results = []
    
    for kid, group in df.groupby('Kid'):
        mean = group[col].mean()
        median = group[col].median()
        q=group[col].quantile(0.75)
        std = group[col].std()
        if opt==1:
            height = lod
        elif opt==2:
            height = max(lod, mean)
        elif opt==3:
            height = max(lod, median)
        else:
            height = max(lod, q)
    
        peaks, properties = sig.find_peaks(group[col], height=height, prominence=0, wlen=5, threshold=0.05)
        peak_rows = group.iloc[peaks]
        
        for i, (idx, row_data) in enumerate(peak_rows.iterrows()):
            
            result = {
                'Kid': kid,
                'height': height,
                'Timepoint': row_data['Timepoint'],   
                col: row_data[col],              
                'PeakHeight': properties['peak_heights'][i],  
                'Prominence': properties['prominences'][i],  
                'log2FoldChange': row_data['log2FoldChange'],
                'LeftBase': properties['left_bases'][i],
                'RightBase': properties['right_bases'][i],
                'LeftThreshold': properties['left_thresholds'][i],
                'RightThreshold': properties['right_thresholds'][i]
            }
            
            if (row_data['log2FoldChange'] < 1 and row_data['log2FoldChange'] == result['Prominence'] and 
                result['RightBase'] - result['LeftBase'] == 2):
                result['peak'] = False
                result['falsetype'] = 1
            elif (row_data['log2FoldChange'] < 1 and result['Prominence'] == result['RightThreshold']):
                result['peak'] = False
                result['falsetype'] = 2
            elif (row_data['log2FoldChange'] < 1):
                result['peak'] = False
                result['falsetype'] = 3
            else:
                result['peak'] = True
                result['falsetype'] =0
            results.append(result)
    
    results = pd.DataFrame(results) 
    return results

def plot_peaks_for_random_kids(data, num_kids=10, random_state=1):

    """
    Plots the original data and detected peaks for a random selection of kids.

    Args:
        data (pd.DataFrame): The original DataFrame containing the data.
        peak_data (pd.DataFrame): The DataFrame containing the detected peaks.
        col (str): The column name to plot (e.g., 'qPCR').
        num_kids (int): The number of random kids to plot.
        random_state (int): The random state for reproducibility.
    """
    # Set seaborn style
    sns.set_theme(style="ticks")

    # Filter random kids
    filter_kid = data['Kid'].drop_duplicates().sample(n=num_kids, random_state=random_state)

    # Filter the peak data for the selected kids
    filtered_data = data[data['Kid'].isin(filter_kid)]

    
    palette = sns.color_palette("tab10", num_kids)

    g= sns.relplot(data=filtered_data, x="Timepoint", y="log2_qPCR",col="Kid", col_wrap=5,kind='line',palette=palette)
    g.savefig('..\\plots\\random_kids.pdf') 

def find_peaks_to(df, col, lod):
    """
    Finds peaks in a DataFrame column using findpeaks.stats.topology() and creates a new DataFrame with the peak information.
    Parameters:
    - df (pandas.DataFrame): The input DataFrame.
    - col (str): The column name in the DataFrame to find peaks in.
    - lod (float): The limit of detection for the peaks.
    - num_std (int, optional): The number of standard deviations to consider when calculating the threshold. Default is 1.
    Returns:
    - pandas.DataFrame: A new DataFrame containing the peak information, with columns 'Kid', 'Timepoint', and the specified column name.
    """
    
    results1 = []
    results2 = []
    

    for kid, group in df.groupby('Kid'):
            if (group is None) or (group.empty) or (len(group) == 1):
                continue
            
            # Use findpeaks to detect peaks
            fp = findpeaks(method='topology', lookahead=1)
            peaks = fp.peaks1d(X=group[col], method='topology')
            #fp.plot_persistence(figsize=(20, 8), fontsize_ax1=14, fontsize_ax2=14, xlabel='x-axis', ylabel='qP')
    

            if (peaks is None) or ('df' not in peaks) or (peaks['df'] is None) or (len(peaks['df']) == 0):
                continue
            
            peaks['df'].rename(columns={'y': col}, inplace=True)
            peaks['df']['Timepoint'] = group['Timepoint'].values
            peaks['df']['Kid'] = kid
            peaks['df']['peak'] = peaks['df']['peak'].astype(bool)

            peaks2 = peaks['df'][(peaks['df'][col] >= lod) & (peaks['df']['peak'] == True) &(peaks['df']['valley'] == False)]
            # Append the DataFrame to the results list
            results1.append(peaks['df'])
            results2.append(peaks2)
         
    final_df1 = pd.concat(results1, ignore_index=True)
    final_df2 = pd.concat(results2, ignore_index=True)

    return final_df1, final_df2

def mergedf(df1, df2):
    merged_df = pd.concat([df1, df2], ignore_index=True)
    df = merged_df[['Kid', 'Method', 'log2_qPCR', 'Timepoint', 'peak', 'falsetype']]		

    return df



def plot_peaks(data):
    """
    Function to generate a percentage stacked bar plot showing the number of peaks by method and timepoint.
    
    Parameters:
    data (DataFrame): DataFrame containing peak data with columns 'Kid', 'Timepoint', 'local', 'topology', 'Method'.
    """
    sns.set_theme(style="ticks")
    grouped_df = data.groupby(['Timepoint', 'Method']).size().reset_index(name='NumberOfPeaks')
    total_peaks_df = grouped_df.groupby('Timepoint')['NumberOfPeaks'].sum().reset_index(name='TotalPeaks')

    merged_grouped_df = pd.merge(grouped_df, total_peaks_df, on='Timepoint')
    merged_grouped_df['Percentage'] = (merged_grouped_df['NumberOfPeaks'] / merged_grouped_df['TotalPeaks']) * 100

    pivot_df = merged_grouped_df.pivot(index='Timepoint', columns='Method', values='Percentage').fillna(0).reset_index()
    melted_df = pivot_df.melt(id_vars='Timepoint', value_vars=['topology', 'local', 'both'], var_name='Method', value_name='Percentage')
    plt.figure(figsize=(14, 14))

    plot = (
        so.Plot(melted_df, x="Timepoint", y="Percentage", color="Method")
        .add(so.Bar(), so.Stack())
        .scale(color="Set2_r")
        .label(x="Timepoint (weeks)", y="Percentage (%)")
        .limit(x=(2, 24), y=(0, 100))  
    )

    fig, ax = plt.subplots()
    plot.on(ax).plot()

    
    ax.set_xticks(range(2, 24, 2))

    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles[2:], loc=10, bbox_to_anchor=(0.5, -0.25), ncol=2, frameon=False, fontsize=14)

    plt.show()




def create_pivot_df(final_df):
    """
    Function to create a pivot DataFrame with 'Method' column and preserve the 'peak' column.
    
    Parameters:
    final_df (DataFrame): DataFrame containing the columns 'Kid', 'Method', 'log2_qPCR', 'Timepoint', and 'peak'.
    
    Returns:
    DataFrame: A DataFrame with 'Kid', 'Timepoint', 'Method', and 'peak' columns.
    """
    pivot_df = final_df.pivot_table(index=['Kid', 'Timepoint'], columns='Method', values=['log2_qPCR', 'peak', 'falsetype'], aggfunc='first')
    pivot_df.columns = ['_'.join(col).strip() for col in pivot_df.columns.values]
    pivot_df = pivot_df.reset_index()
    pivot_df['Method'] = pivot_df.apply(
    lambda row: 'both' if pd.notna(row['log2_qPCR_local']) and pd.notna(row['log2_qPCR_topology']) else 
                ('topology' if pd.notna(row['log2_qPCR_topology']) else 'local'), axis=1)
    pivot_df['log2_qPCR'] = pivot_df[['log2_qPCR_local', 'log2_qPCR_topology']].bfill(axis=1).iloc[:, 0]
    pivot_df['peak'] = pivot_df[['peak_local', 'peak_topology']].bfill(axis=1).iloc[:, 0]
    pivot_df=pivot_df.rename(columns={'falsetype_local':'falsetype'}).drop(columns=[
    'log2_qPCR_local', 
    'log2_qPCR_topology', 
    'peak_local', 
    'peak_topology'])
    pivot_df['falsetype']=pivot_df['falsetype'].fillna(0)

    return pivot_df


def plot_heatmap(final_df):
    """
    Function to generate a heatmap showing the Method values for each Kid and Timepoint.
    
    Parameters:
    final_df (DataFrame): DataFrame containing the columns 'Kid', 'Method', 'log2_qPCR', and 'Timepoint'.
    """
    pivot_df = create_pivot_df(final_df)
    Method_mapping = {'topology': 0, 'local': 1, 'both': 2}
    pivot_df['Method_num'] = pivot_df['Method'].map(Method_mapping)
    heatmap_data = pivot_df.pivot(index='Kid', columns='Timepoint', values='Method_num')
    cmap = sns.color_palette(['#fcdc4c','#de79f2','#f55953'], as_cmap=True)
    plt.figure(figsize=(12, 8))
    ax = sns.heatmap(heatmap_data, cmap=cmap, linewidths=.2, linecolor='gray', cbar=False, annot=False, fmt='')

    plt.ylabel('Kids', fontsize=16)
    plt.xlabel('Timepoint (weeks)', fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(rotation=360, fontsize=12)
    plt.grid(False)

    plt.tight_layout()
    plt.show()
    
def plot_levels(df, col):
    heatmap_data = df.pivot(index='Kid', columns='Timepoint', values=col)
    
    cmap = sns.color_palette("mako", as_cmap=True).reversed()
    plt.figure(figsize=(12, 8))
    ax = sns.heatmap(
        heatmap_data, cmap=cmap, linewidths=.05, linecolor='gray',  # Ajustar el grosor de las líneas de la cuadrícula
        cbar_kws={'label': col, 'shrink': 0.5}  # Ajustar el tamaño de la barra de color
    )
    
    plt.ylabel('Kids', fontsize=20)
    plt.xlabel('Timepoint (weeks)', fontsize=20)
    plt.xticks(fontsize=12)
    plt.yticks(rotation=360, fontsize=12)
    plt.grid(False)
    
    plt.tight_layout()
    plt.show()




    

