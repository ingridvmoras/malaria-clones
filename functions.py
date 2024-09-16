
from matplotlib import pyplot as plt
import numpy as np
import scipy.signal as sig
import pandas as pd
import seaborn as sns
from findpeaks import findpeaks
from matplotlib.colors import ListedColormap
from matplotlib.backends.backend_pdf import PdfPages

def rollavg_convolve(a, n):
    'scipy.convolve'
    assert n % 2 == 1
    if len(a) < n:
        return np.full_like(a, np.nan)  # Return NaNs if the array is too short
    convolved = sig.convolve(a, np.ones(n, dtype='float') / n, 'same')
    pad_size = n // 2
    return np.pad(convolved[pad_size:-pad_size], (pad_size, pad_size), mode='edge')


# def find_peaks_dx(df, col):
#     """
#     Detects peaks in parasitemia data based on derivatives.

#     Args:
#         df: A pandas DataFrame containing the parasitemia data.
#         col: The column name to detect peaks in.

#     Returns:
#         A pandas DataFrame containing the peak information for each child.
#     """

#     # Group by child
#     grouped_data = df.groupby('Kid')

#     # Calculate derivatives and identify peaks for each child
#     results = []
#     for kid, group in grouped_data:
#         if len(group) < 2:
#             continue  # Skip groups with less than 2 elements

#         # Calculate first derivative
#         first_derivative = np.gradient(group[col])

#         # Calculate second derivative
#         second_derivative = np.gradient(first_derivative)

#         # Identify peaks where first derivative changes sign and second derivative is negative
#         peaks = np.where(np.diff(np.sign(first_derivative)) == -2)[0] + 1
#         peaks = peaks[second_derivative[peaks] < 0]

#         if len(peaks) == 0:
#             continue  # Skip if no peaks are found

#         # Extract rows corresponding to the peak indices
#         peak_rows = group.iloc[peaks]

#         # Append the results with Kid, Timepoint, and qPCR
#         for _, row in peak_rows.iterrows():
#             results.append({
#                 'Kid': kid,
#                 'Timepoint': row['Timepoint'],
#                 col: row[col]
#             })

#     return pd.DataFrame(results)


def find_peaks_lm(df, col, lod, num_std=1):
    """
    Finds peaks in a DataFrame column and creates a new DataFrame with the peak information.
    Parameters:
    - df (pandas.DataFrame): The input DataFrame.
    - col (str): The column name in the DataFrame to find peaks in.
    - lod (float): The limit of detection for the peaks.
    - num_std (int, optional): The number of standard deviations to consider when calculating the threshold. Default is 1.
    Returns:
    - pandas.DataFrame: A new DataFrame containing the peak information, with columns 'Kid', 'Timepoint', and the specified column name.
    """
    
    results = []
    for kid, group in df.groupby('Kid'):
        # Calcular el umbral, considerando el LOD y la desviación estándar
        mean = group[col].mean()
        std = group[col].std()
        threshold = max(lod, mean + num_std * std)

        # Encontrar los picos
        peaks, _ = sig.find_peaks(group[col], height=threshold)
        
        # Extract rows corresponding to the peak indices
        peak_rows = group.iloc[peaks]
        
        # Append the results with Kid, Timepoint, and qPCR
        for _, row in peak_rows.iterrows():
            results.append({
                'Kid': kid,
                'Timepoint': row['Timepoint'],
                col: row[col],
                'Threshold': threshold
            })

    return pd.DataFrame(results)



def plot_peaks_for_random_kids(data, peak_data, col, num_kids=10, random_state=1):
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
    sns.set_theme(style="whitegrid")

    # Filter random kids
    filter_kid = data['Kid'].drop_duplicates().sample(n=num_kids, random_state=random_state)

    # Filter the peak data for the selected kids
    filtered_peak = peak_data[peak_data['Kid'].isin(filter_kid)]
    filtered_data = data[data['Kid'].isin(filter_kid)]

    
    palette = sns.color_palette("tab10", num_kids)

    # Plot the data
    plt.figure(figsize=(12, 6))

    # Plot the detected peaks
    plt.scatter(filtered_peak['Timepoint'], filtered_peak[col], color='red', label='Detected Peaks', zorder=5)

    # Plot the original data for each kid with different colors
    for i, kid in enumerate(filter_kid):
        kid_data = filtered_data[filtered_data['Kid'] == kid]
        plt.plot(kid_data['Timepoint'], kid_data[col], marker='o', linestyle='-', color=palette[i], label=f'Kid {kid}', zorder=1)

    plt.xlabel('Timepoint (weeks)')
    plt.ylabel(f'Parasitemia (Log-transformed {col})')
    plt.legend(title='Child', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.title(f'Detected Peaks in {col} Data for Selected Kids')

    plt.xticks(ticks=filtered_data['Timepoint'].unique())
    plt.tight_layout()

    plt.show()

# Example usage
# plot_peaks_for_random_kids(first_qpcr, peak_data_dxqpcr, 'qPCR')


def find_peaks_to(df, col, lod, num_std=1):
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
    
    #with PdfPages('peaks_plots.pdf') as pdf:
    for kid, group in df.groupby('Kid'):
            if (group is None) or (group.empty) or (len(group) == 1):
                continue
            
            # Use findpeaks to detect peaks
            fp = findpeaks(method='topology', lookahead=1)
            peaks = fp.peaks1d(X=group[col], method='topology')
               
            # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
            # ax1, ax2 = fp.plot_persistence(figsize=(20, 8), fontsize_ax1=14, fontsize_ax2=14, xlabel='x-axis', ylabel='y-axis')
            # fig.suptitle(f'Persistence Plot for Kid {kid}')
            # pdf.savefig(fig)
            # plt.close(fig)

            if (peaks is None) or ('df' not in peaks) or (peaks['df'] is None) or (len(peaks['df']) == 0):
                continue
            
            peaks['df'].rename(columns={'y': col}, inplace=True)
            peaks['df']['Timepoint'] = group['Timepoint'].values
            peaks['df']['Kid'] = kid
            
            # Ensure 'peak' column is boolean
            peaks['df']['peak'] = peaks['df']['peak'].astype(bool)
            
            # Filter the rows where 'qPCR' is above the threshold or 'peak' is True
            peaks2 = peaks['df'][(peaks['df'][col] >= 2) & (peaks['df']['peak'] == True)]
            
            # Append the DataFrame to the results list
            results1.append(peaks['df'])
            results2.append(peaks2)
         
    
    # Concatenate all DataFrames in the results list into a single DataFrame
    final_df1 = pd.concat(results1, ignore_index=True)
    final_df2 = pd.concat(results2, ignore_index=True)

    return final_df1, final_df2

def plot_peaks(peaks, peaks_lm_qpcr):
    """
    Function to generate a bar plot showing the number of peaks by method and timepoint.
    
    Parameters:
    peaksto1 (DataFrame): DataFrame containing peak data for the 'topology' method.
    peaks_lm_qpcr (DataFrame): DataFrame containing peak data for the 'local' method.
    """
    sns.set_theme(style="ticks")
    merged_df = pd.concat([peaks, peaks_lm_qpcr], ignore_index=True)

    final_df = merged_df[['Kid', 'Method', 'qPCR', 'Timepoint']]
    print(final_df)

    grouped_df = final_df.groupby(['Method', 'Timepoint']).size().reset_index(name='NumberOfPeaks')
    timepoint_order = sorted(grouped_df['Timepoint'].unique())  # Sort the Timepoints
    method_order = ['topology', 'local']  
    
    
    plt.figure(figsize=(12, 8))
    ax = sns.barplot(data=grouped_df, x='Timepoint', y='NumberOfPeaks', hue='Method', palette= 'Paired', order=timepoint_order, hue_order=method_order)
    ax.set_xlabel('Timepoint (weeks)', fontsize=14)
    ax.set_ylabel('Number of Peaks', fontsize=14)
    plt.xticks(rotation=90, fontsize=12)
    plt.yticks(fontsize=12)
    sns.despine(offset=10, trim=True)
    plt.tight_layout()
    plt.show()
    
    return final_df



def create_pivot_df(final_df):
    """
    Function to create a pivot DataFrame with 'identify_by' column.
    
    Parameters:
    final_df (DataFrame): DataFrame containing the columns 'Kid', 'Method', 'qPCR', and 'Timepoint'.
    
    Returns:
    DataFrame: A DataFrame with 'Kid', 'Timepoint', and 'identify_by' columns.
    """
    pivot_df = final_df.pivot_table(index=['Kid', 'Timepoint'], columns='Method', values='qPCR', aggfunc='first')
    pivot_df['identify_by'] = pivot_df.apply(lambda row: 'both' if pd.notna(row['topology']) and pd.notna(row['local']) else ('topology' if pd.notna(row['topology']) else 'local'), axis=1)
    pivot_df = pivot_df.reset_index()
    
    return pivot_df



def plot_heatmap(final_df):
    """
    Function to generate a heatmap showing the identify_by values for each Kid and Timepoint.
    
    Parameters:
    final_df (DataFrame): DataFrame containing the columns 'Kid', 'Method', 'qPCR', and 'Timepoint'.
    """
    pivot_df = create_pivot_df(final_df)
    identify_by_mapping = {'topology': 0, 'local': 1, 'both': 2}
    pivot_df['identify_by_num'] = pivot_df['identify_by'].map(identify_by_mapping)
    heatmap_data = pivot_df.pivot(index='Kid', columns='Timepoint', values='identify_by_num')
    cmap = sns.color_palette("Paired", as_cmap=True)
    plt.figure(figsize=(20, 18))
    ax = sns.heatmap(heatmap_data, cmap=cmap, linewidths=.2, linecolor='gray', cbar=True, annot=False, fmt='')

    plt.ylabel('Kids', fontsize=16)
    plt.xlabel('Timepoint (weeks)', fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(rotation=360, fontsize=12)
    plt.grid(False)

    plt.tight_layout()
    plt.show()
    
def plot_levels(df,col):
    heatmap_data = df.pivot(index='Kid', columns='Timepoint', values= col)
    
    cmap = sns.color_palette("ch:start=.2,rot=-.3", as_cmap=True)
    plt.figure(figsize=(20, 18))
    ax = sns.heatmap(heatmap_data, cmap=cmap, linewidths=.2, linecolor='gray', cbar_kws={'label': col })
    
    plt.ylabel('Kids', fontsize=16)
    plt.xlabel('Timepoint (weeks)', fontsize=16)
    plt.xticks(fontsize=12)
    plt.yticks(rotation=360, fontsize=12)
    plt.grid(False)
    
    plt.tight_layout()
    plt.show() 