import pandas as pd
import numpy as np

import pandas as pd
from datetime import datetime

def chunk_analysis(df_sorted,df, max_zeroes):
    """
    timepoint chunk analysis for each haplotype in the dataframe.

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





# Load data
df = pd.read_csv('..\\outcome\\peaks.csv')
df['cluster_name'] = df['cluster_name'].str.replace("_", "")

# Prepare the data
original_unique = df.copy()  # To add metadata to the length file

df['count'] = 1  

df['Sample'] = df['Kid'].astype(str) + "_" + df['cluster_name']
df = df.drop_duplicates(subset=['Sample', 'Timepoint', 'count'])
df = df.pivot(index='Timepoint', columns='Sample', values='count').fillna(0).astype(int)

df_sorted = df.sort_index().reset_index()

# Define the list of max_zeroes values
max_zeroes_list = [0, 1, 2, 3]

# Get the current date
current_date = datetime.now().strftime("%Y-%m-%d")

# Loop through max_zeroes_list and save results
for max_zeroes in max_zeroes_list:
    result_df = chunk_analysis(df_sorted,original_unique, max_zeroes)
    file_name = f"..\\outcome\\infection_duration_analysis_{max_zeroes}_zeroes_{current_date}.csv"
    result_df.to_csv(file_name, index=False)
