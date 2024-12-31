
from matplotlib import pyplot as plt
import numpy as np
import scipy as sci
import pandas as pd
import seaborn as sns
import seaborn.objects as so
from matplotlib.colors import ListedColormap
from matplotlib.backends.backend_pdf import PdfPages



# def rollavg_convolve(a, n):
#     'scipy.convolve'
#     assert n % 2 == 1
#     if len(a) < n:
#         return np.full_like(a, np.nan) 
#     convolved = sig.convolve(a, np.ones(n, dtype='float') / n, 'same')
#     pad_size = n // 2
#     return np.pad(convolved[pad_size:-pad_size], (pad_size, pad_size), mode='edge')

def mergedf(dfs:list,col:str, id:str):
    """_summary_

    Args:
        dfs (list): List of dataframes 
        col (str): Column name to merge with parasitemia data, could be 'log2_qPCR', 'qPCR', etc.
        id (str): Column name to merge with the id of individuals, could be 'Kid', etc.

    Returns:
        _type_: pd.Dataframe with the merged data 
    """    
    merged = pd.concat(dfs, ignore_index=True)
    df = merged.filter(items=[id, 'Timepoint', col, 'Method', 'peak', 'valley'])

    return df

def sample_random(df, id, num_samples=1, random_state=41):
    """
    Groups the DataFrame by 'id' and obtains a random sample from each group.

    Args:
        df (pd.DataFrame): The original DataFrame containing the data.
        num_samples (int): The number of random samples to obtain for each 'id'.
        random_state (int): The random state for reproducibility.

    Returns:
        pd.DataFrame: The DataFrame with the random sample.
    """
  
    sampled_ids = df[id].drop_duplicates().sample(n=num_samples, random_state=random_state)
    sample =  df[df[id].isin(sampled_ids)]
    return sample

def merge_peak_data(data: pd.DataFrame, df_peak: pd.DataFrame,method, id='Kid') -> pd.DataFrame:
    """
    Converts the 'peak' column of df_peak to 1 if True and 0 if False,
    and merges df_peak into the DataFrame data using the id and Timepoint columns as join keys.

    Args:
        data (pd.DataFrame): Original DataFrame.
        df_peak (pd.DataFrame): DataFrame with the 'peak' column.

    Returns:
        pd.DataFrame: Combined DataFrame with the 'peak' column merged.
    """
    df_peak['peak'] = df_peak['peak'].astype(int)
    if method == 'LM' or method == 'TPH':
        df_peak['valley'] = df_peak['valley'].astype(int)
        merged_data = pd.merge(data, df_peak[[id, 'Timepoint', 'peak','valley']], on=[id, 'Timepoint'], how='left')
    else:
        merged_data = pd.merge(data, df_peak[[id, 'Timepoint', 'peak']], on=[id, 'Timepoint'], how='left')

    return merged_data
