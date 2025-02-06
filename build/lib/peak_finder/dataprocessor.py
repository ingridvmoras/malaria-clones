# data_processor.py
import pandas as pd
import numpy as np

class DataProcessor:
    def __init__(self, file_path):
        self.file_path = file_path

    def load_data(self, filepath):
        return pd.read_csv(filepath)

    def preprocess_data(self, data, col='qPCR',id='Kid'):
        may_timepoint = 13
        data = data.replace('May13', may_timepoint, regex=True).replace('May12', 0, regex=True)
        data['Timepoint'] = pd.to_numeric(data['Timepoint'])

        # Filter out May12 and May13 timepoints
        filtered_data = data[~data['Timepoint'].isin([may_timepoint, 0])].copy()
        first_data = filtered_data.drop_duplicates(subset=['Timepoint', id], keep='first').copy()
        first_data = first_data.sort_values(by=[id, 'Timepoint'])

        # NA values are filled with the mean of the qPCR values for each kid. qPCR values lower than 1 are set to 1.
        first_data[col] = first_data[col].fillna(first_data.groupby(id)[col].transform('mean'))
        first_data[col] = first_data[col].clip(lower=1)

        return first_data, data

    def process_data(self, data, col='qPCR', id='Kid'):
       
        # Log transform qPCR counts for more normal distribution
        data[f'log2_{col}'] = np.log2(data[col])
        data[f'log10_{col}'] = np.log10(data[col])
        data['FoldChange'] = data.groupby(id)[f'{col}'].diff()
        data['log2FoldChange'] = data.groupby(id)[f'log2_{col}'].diff()
        data['log10FoldChange'] = data.groupby(id)[f'log10_{col}'].diff()
        data = data[data[id].isin(data[id].value_counts()[data[id].value_counts() >= 3].index)]
        data = data.sort_values(by=[id, 'Timepoint'])
        num_ids = data[id].nunique()
        print(f"Number of parasitaemia time series for peak detection: {num_ids}")
        return data, num_ids
    
    


