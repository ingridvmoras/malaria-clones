
from matplotlib import pyplot as plt
import numpy as np
import scipy.signal as sig
import pandas as pd
import seaborn as sns

from detecta import detect_peaks

def rollavg_convolve(a, n):
    'scipy.convolve'
    assert n % 2 == 1
    if len(a) < n:
        return np.full_like(a, np.nan)  # Return NaNs if the array is too short
    convolved = sig.convolve(a, np.ones(n, dtype='float') / n, 'same')
    pad_size = n // 2
    return np.pad(convolved[pad_size:-pad_size], (pad_size, pad_size), mode='edge')


def find_peaks_dx(df, col):
    """
    Detects peaks in parasitemia data based on derivatives.

    Args:
        df: A pandas DataFrame containing the parasitemia data.
        col: The column name to detect peaks in.

    Returns:
        A pandas DataFrame containing the peak information for each child.
    """

    # Group by child
    grouped_data = df.groupby('Kid')

    # Calculate derivatives and identify peaks for each child
    results = []
    for kid, group in grouped_data:
        if len(group) < 2:
            continue  # Skip groups with less than 2 elements

        # Calculate first derivative
        first_derivative = np.gradient(group[col])

        # Calculate second derivative
        second_derivative = np.gradient(first_derivative)

        # Identify peaks where first derivative changes sign and second derivative is negative
        peaks = np.where(np.diff(np.sign(first_derivative)) == -2)[0] + 1
        peaks = peaks[second_derivative[peaks] < 0]

        if len(peaks) == 0:
            continue  # Skip if no peaks are found

        # Extract rows corresponding to the peak indices
        peak_rows = group.iloc[peaks]

        # Append the results with Kid, Timepoint, and qPCR
        for _, row in peak_rows.iterrows():
            results.append({
                'Kid': kid,
                'Timepoint': row['Timepoint'],
                col: row[col]
            })

    return pd.DataFrame(results)


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



def plot_peaks_for_random_kids(data, peak_data, col, num_kids=10, random_state=42):
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


