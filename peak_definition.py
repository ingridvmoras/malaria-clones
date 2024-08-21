import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sci
import statannot as sta  # import add_stat_annotation

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

# Select the first qPCR value for each timepoint per kid (it does not matter which one, 
# since parasitemia is the same for each value in the same timepoint)

first_qpcr = filtered_dat.groupby(['Timepoint', 'Kid']).first().reset_index()
# Randomly choose 10 unique kids
random_kids = first_qpcr['Kid'].drop_duplicates().sample(n=10, random_state=42)

# Filter the data for the randomly chosen kids
random_kids_data = first_qpcr[first_qpcr['Kid'].isin(random_kids)]

# Plot the histogram for the randomly chosen kids
plt.figure(figsize=(12, 6))
for kid in random_kids:
    kid_data = random_kids_data[random_kids_data['Kid'] == kid]
    plt.plot(kid_data['Timepoint'], kid_data['qPCR'], marker='o', linestyle='-', label=f'Kid {kid}')

plt.xlabel('Timepoint (weeks)')
plt.ylabel('Parasitemia (Log-transformed qPCR)')
plt.legend(title='Child', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()