
import random
import numpy 
import pandas as pd
from datetime import datetime
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

def sample_random(df, id, num_samples=1, random_state=42):
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

def chunk_analysis(df_sorted,df, max_zeroes):
    """
    Timepoint chunk analysis for each haplotype in the dataframe.

    """
    df_length = []

    # Loop through each haplotype and find chunks of contiguous timepoints
    for name in df_sorted.columns[1:]:
        values = df_sorted[name].values
        start = None
        end = None
        num_zeroes = 0
        chunks = []
        last_one_seen = None

        for i, value in enumerate(values):
            if value == 0:
                if start is not None:
                    num_zeroes += 1
                    if num_zeroes > max_zeroes:
                        end = i - 1 - max_zeroes
                        chunks.append((start, end))
                        start = None
                        end = None
                        num_zeroes = 0
            else:
                if start is None:
                    start = i
                num_zeroes = 0
                last_one_seen = i

        # Handle the last value in the column
        if start is not None:
            if values[-1] == 0:
                end = last_one_seen
            else:
                end = len(values) - 1
            chunks.append((start, end))
        elif values[-1] == 1:
            chunks.append((len(values) - 1, len(values) - 1))

        # Extract Kid and cluster_name from haplotype name
        Kid, cluster_name = name.split("_")
        for chunk in chunks:
            start = df_sorted['Timepoint'].iloc[chunk[0]]
            end = df_sorted['Timepoint'].iloc[chunk[1]]
            df_length.append([Kid, cluster_name, start, end])

    # Transform the list of chunks into a DataFrame
    df_length = pd.DataFrame(df_length, columns=["Kid", "cluster_name", "Start", "End"])
    merged = pd.merge(df,df_length, how='left', on=['Kid', 'cluster_name'])
    merged['first_appearance'] = (merged['Start'] == merged['Timepoint']).astype(int)
    merged['first_TP'] = (merged['Start'] == merged['Timepoint']).astype(bool)
    return merged




def first_appearances(df: pd.DataFrame, id: str, method: str, output_dir:str) -> pd.DataFrame:
    """
    Analyzes the first appearances of clusters in a DataFrame and saves the results to CSV files.
    This function processes a DataFrame to identify the first appearances of clusters based on a given ID and method.
    It generates multiple CSV files with the analysis results, each corresponding to different maximum zero values.
    Args:
        df (pd.DataFrame): The input DataFrame containing the data to be analyzed.
        id (str): The column name in the DataFrame that contains the unique identifier for each sample.
        method (str): The method used for the analysis (LM,TPH,S1), which will be included in the output file names.
        output_dir (str): The directory where the output CSV files will be saved.
    Returns:
        pd.DataFrame: The processed DataFrame with the first appearances of clusters.
    
    """
    df['cluster_name'] = df['cluster_name'].str.replace("_", "")
    original_unique = df.copy()  # To add metadata to the length file
    df['count'] = 1  
    df['Sample'] = df[id].astype(str) + "_" + df['cluster_name']
    df = df.drop_duplicates(subset=['Sample', 'Timepoint', 'count'])
    df = df.pivot(index='Timepoint', columns='Sample', values='count').fillna(0).astype(int)
    df_sorted = df.sort_index().reset_index()
    
    # Define the list of max_zeroes values
    max_zeroes_list = [0, 1, 2, 3]
    
    # Get the current date
    current_date = datetime.now().strftime("%Y-%m-%d")
    
    # Loop through max_zeroes_list and save results
    for max_zeroes in max_zeroes_list:
        result_df = chunk_analysis(df_sorted, original_unique, max_zeroes)
        file_name = f"{output_dir}/infection_duration_analysis_{max_zeroes}_zeroes_{current_date}_{method}.csv"
        result_df.to_csv(file_name, index=False)

 
    
    
def run_simulation(df,id='Kid',col='qPCR', rounds=100, weighted=True):
    """
    This function simulates the first appearance of clones in a dataset based on the given parameters.
    Args:
        df (pd.DataFrame): The input DataFrame containing the data to be analyzed.
        id (str): The column name in the DataFrame that contains the unique identifier for each sample.
        col (str): The column name in the DataFrame that contains the data to be analyzed.
        rounds (int): The number of simulation rounds to run.
        weighted (bool): A flag indicating whether to use weighted sampling.
    """
    #Code modified from apperances_Sept2023.R by Manuela Carrasquilla
    
    # keep track of which clones are seen at each time point
    # includes the previous timepoint
    # also keeping track of the frequency of each clone
    c_tp = {}
    for tp in df['Timepoint'].unique():
        v =df[df['Timepoint'].isin([tp, tp-1])]['cluster_name'].value_counts() / df[df['Timepoint'].isin([tp, tp-1])].shape[0]
        c_tp[tp] = (list(v.index), v.values)
    
    # run simulations
    # do it cycling through each timepoint first
    # then kid
    # then each simulation
    # this is probably the fastest way to do it
    res = []
    for tp in sorted(df['Timepoint'].unique()):
        print(f'simulating timepoint {tp}')
        
        # select time point from original dataframe
        t = df[df['Timepoint'] == tp]

        # cycle through each id in this time point
        for ind in t[id].unique():
            # select id from this timepoint
            tt = t[t[id] == ind]
            
            # keep all the values that are not gonna be randomized
            # in a variable for convenience
            values = tt.iloc[0]
            
            # run the randomizations
            for i in range(rounds):
                # pick random clones
                while True:
                    # weighted: the proportion of clones is taken into account
                    if weighted:
                        clones = random.choices(c_tp[tp][0], k=values['COI'],
                                        weights=c_tp[tp][1])
                    # unweighted: each clone is equally probable
                    else:
                        clones = random.choices(c_tp[tp][0], k=values['COI'])
                    
                    # avoid picking the same clone multiple times
                    # inefficient but whatever
                    if len(set(clones)) == values['COI']:
                        break
                
                # save the data in a list
                for clone in clones:
                    res.append((ind, tp, clone,
                                values['COI'],
                                values[col],
                                values['peak'],
                                i))

    # create the final dataframe
    r = pd.DataFrame(res,
                     columns=[id, 'Timepoint', 'cluster_name',
                              'COI', col,
                              'peak', 'simulation'])
    
    # keep track of when each clone in each kid and simulation was first observed
    r['first_appearance'] = 0
    r = r.set_index([id, 'cluster_name', 'simulation', 'Timepoint']).sort_index()
    r.loc[r.reset_index().groupby([id, 'cluster_name', 'simulation'])['Timepoint'].min().reset_index().values,
          'first_appearance'] = 1
    r = r.reset_index()
    r = r.sort_values(['Timepoint', id, 'simulation', 'cluster_name'])
    
    return r

