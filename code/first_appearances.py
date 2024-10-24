import pandas as pd
import numpy as np

# Load raw data from the CSV file
df = pd.read_csv('..\\outcome\\peaks.csv')

# Format and clean the data
df = df.sort_values(by='Date')
original_df = df.copy()
df['count'] = 1  # Add counter

original_unique = df.copy()  # To add metadata to the length file

df['cluster_name'] = df['cluster_name'].str.replace("_", "")

df['Sample'] = df['Kid'].astype(str) + "_" + df['cluster_name']
df['Sample'] = df['Sample'].str.replace("\\.", "", regex=True)
df = df.drop_duplicates(subset=['Sample', 'Timepoint', 'count'])
df = df.pivot(index='Timepoint', columns='Sample', values='count').fillna(0).astype(int)

df_sorted = df.sort_index().reset_index()

# Timepoint chunk analysis
max_zeroes = 2  # Adjust as needed for the desired number of skips between positive timepoints
df_length = []

# Loop through each haplotype and find chunks of contiguous timepoints
for name in df_sorted.columns[1:]:
    values = df_sorted[name].values
    start = None
    end = None
    num_zeroes = 0
    chunks = []
    last_one_seen = None

    for i in range(len(values)):
        if i != len(values) - 1:
            value = values[i]
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

    value = values[-1]
    if start is not None:
        if value == 0:
            end = last_one_seen
            chunks.append((start, end))
        else:
            end = len(values) - 1
            chunks.append((start, end))
    elif value == 1:
        chunks.append((len(values) - 1, len(values) - 1))

    # Extract Kid and cluster_name from haplotype name
    Kid, cluster_name = name.split("_")
    for chunk in chunks:
        start = df_sorted['Timepoint'].iloc[chunk[0]]
        end = df_sorted['Timepoint'].iloc[chunk[1]]
        df_length.append([Kid, cluster_name, start, end])

# Transform the list of chunks into a DataFrame
df_length = pd.DataFrame(df_length, columns=["Kid", "Haplotype", "Start", "End"])

# Map original dates to the identified chunks
df_length['Start'] = df_length['Start'].replace({44: "May13", 0: "May12"})
df_length['End'] = df_length['End'].replace({44: "May13", 0: "May12"})

# DataFrame to a CSV file
df_length.to_csv('..\\outcome\\infection_duration_analysis.csv', index=False)