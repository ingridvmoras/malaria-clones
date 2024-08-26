import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sci
import scipy.stats as stats
import statannot as sta 

sns.set_style("white")
sns.set_context("talk")

# Load the data
dat = pd.read_csv('dataset_Kalifabougou.csv', index_col=0)

# Replace May12 and May13 timepoints with may_timepoint
may_timepoint = 22
dat = dat.replace('May13', may_timepoint, regex=True).replace('May12', may_timepoint, regex=True)

# Convert Timepoint to numeric units of 1 week
t_units = 2
dat['Timepoint'] = pd.to_numeric(dat['Timepoint']) * t_units

# Filter out May12 and May13 timepoints
filtered_dat = dat[~dat['Timepoint'].isin([may_timepoint * t_units])].copy()

# Log transform qPCR counts for more normal distribution
filtered_dat.loc[:, 'qPCR'] = np.log(filtered_dat['qPCR'])


#FOLD_CHANGE 

first_qpcr = filtered_dat.dropna().groupby(['Timepoint', 'Kid']).first().reset_index()

# Calculate fold change for each kid at each timepoint
first_qpcr['FoldChange'] = first_qpcr.groupby('Kid')['qPCR'].pct_change(fill_method=None)+1
first_qpcr['FoldChange'] = first_qpcr['FoldChange'].abs()

# Exclude the first timepoint for each kid (to avoid division by zero)
first_qpcr= first_qpcr.groupby('Kid').apply(lambda x: x.iloc[1:],include_groups=True ).reset_index(drop=True)
