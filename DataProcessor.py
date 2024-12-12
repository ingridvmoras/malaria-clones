# data_processor.py
import pandas as pd
import numpy as np

class DataProcessor:
    def __init__(self, file_path):
        self.file_path = file_path

    def load_data(self, filepath):
        return pd.read_csv(filepath)

    # def preprocess_data(self, data):
    #     may_timepoint = 13
    #     data = data.replace('May13', may_timepoint, regex=True).replace('May12', 0, regex=True)
    #     data['Timepoint'] = pd.to_numeric(data['Timepoint'])
    #     filtered_data = data[~data['Timepoint'].isin([may_timepoint, 0])].copy()
    #     first_data = filtered_data.groupby(['Timepoint', 'Kid']).first().reset_index()
    #     first_data['qPCR'] = first_data['qPCR'].fillna(first_data.groupby('Kid')['qPCR'].transform('mean'))
    #     first_data['qPCR'] = first_data['qPCR'].clip(lower=1)
    #     first_data['log2_qPCR'] = np.log2(first_data['qPCR'])
    #     first_data['log10_qPCR'] = np.log10(first_data['qPCR'])
    #     first_data['log2FoldChange'] = first_data.groupby('Kid')['log2_qPCR'].diff()
    #     first_data['log10FoldChange'] = first_data.groupby('Kid')['log10_qPCR'].diff()
    #     return first_data

    def load_and_process_data(self):
        # Load the data
        data = pd.read_csv(self.file_path)
        
        col= 'qPCR' #Replace with the name of the column that contains the parasitemia values
        id= 'Kid' #Replace with the name of the column that contains the individual IDs
        
        # Log transform qPCR counts for more normal distribution
        data[f'log2_{col}'] = np.log2(data[col])
        data[f'log10_{col}'] = np.log10(data[col])
        data['log2FoldChange'] = data.groupby(id)[f'log2_{col}'].diff()
        data['log10FoldChange'] = data.groupby(id)[f'log10_{col}'].diff()

        # Calculate the moving average of Parasitemia for each kid
        data['MovingAverage'] = data.groupby(id)[col].transform(lambda x: x.rolling(window=3, min_periods=1).mean())

        # Exclude the first timepoint for each kid (to avoid division by zero)
        data_ft = data.copy()  # Dataframe with the first timepoint for each kid
        data = data.groupby(id).apply(lambda x: x.iloc[1:]).reset_index(drop=True)
        data['FoldChange'] = data['log2FoldChange'].apply(lambda x: pow(2, x) if pd.notnull(x) else np.nan)

        return data, data_ft

# main.py
#from data_processor import DataProcessor

# Initialize the DataProcessor with the file path
#processor = DataProcessor('..\\data\\dataset_Kalifabougou.csv')

# Load and process the data
#data, data_ft = processor.load_and_process_data()

# Continue with further analysis or visualization