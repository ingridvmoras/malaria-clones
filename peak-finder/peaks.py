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
import utils as f
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
first_data = first_data.sort_values(by=['Kid', 'Timepoint'])

# NA values are filled with the mean of the qPCR values for each kid. qPCR values lower than 1 are set to 1.

first_data['qPCR'] = first_data['qPCR'].fillna(first_data.groupby('Kid')['qPCR'].transform('mean'))
first_data['qPCR'] = first_data['qPCR'].clip(lower=1)

# Log transform qPCR counts for more normal distribution
first_data['log2_qPCR'] = np.log2(first_data['qPCR'])


first_data['log10_qPCR'] = np.log10(first_data['qPCR'])
first_data['log2FoldChange'] = first_data.groupby('Kid')['log2_qPCR'].diff()
first_data['log10FoldChange'] = first_data.groupby('Kid')['log10_qPCR'].diff()
first_data = first_data[first_data['Kid'].isin(first_data['Kid'].value_counts()[first_data['Kid'].value_counts() >= 3].index)]

#Calculate the moving avarage of Parasitemia for each kid
#first_data['MovingAverage'] = first_data.groupby('Kid')['qPCR'].transform(lambda x: f.rollavg_convolve(x, 3))

# Exclude the first timepoint for each kid (to avoid division by zero)
first_data_ft= first_data.copy() #Dataframe with the first timepoint for each kid
first_data= first_data.groupby('Kid').apply(lambda x: x.iloc[1:],include_groups=True ).reset_index(drop=True)
first_data['FoldChange'] = first_data['log2FoldChange'].apply(lambda x: pow(2, x))




