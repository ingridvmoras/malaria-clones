import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sci
import scipy.stats as stats
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

first_qpcr = filtered_dat.groupby(['Timepoint', 'Kid']).first().reset_index().dropna()

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

## First trial  (Threshold-based change rate peak detection)
threshold = 2  # Parasitemia threshold

data = kid_data.copy()

# Calculate the rate of change using first differences
data['rate_of_change'] = data['qPCR'].diff()

# Calculate a moving average (adjust window size as needed)
window_size = 2
data['moving_average'] = data['qPCR'].rolling(window=window_size).mean()

# Identify peaks based on threshold, rate of change, moving average, and trend
peaks = data[(data['qPCR'] > threshold)  & (data['rate_of_change'] < 0) & (data['qPCR'] > data['moving_average'])]

# Group data by peak and non-peak periods
non_peak_data = data[~data.index.isin(peaks.index)]

# # Perform a Mann-Whitney U test
# u_statistic, p_value = stats.mannwhitneyu(peaks['qPCR'], non_peak_data['qPCR'])

# # Print results
# print("U-statistic:", u_statistic)
# print("P-value:", p_value)

# # Calculate Cohen's d
# cohen_d = u_statistic / np.sqrt((len(peaks) + len(non_peak_data)) / 2)

# print("Cohen's d:", cohen_d)

# With this method I can detect peaks based on a threshold, rate of change, moving average, and trend. However, 
# the results are not as clear as I would like. I will try another method to detect peaks :( )

## Second trial (Peak detection using scipy.signal.find_peaks)

peaks, _ = sci.signal.find_peaks(data['qPCR'], height=threshold)
 
plt.figure(figsize=(12, 6))
plt.plot(data['Timepoint'], data['qPCR'], marker='o', linestyle='-', label='Data')
plt.plot(data['Timepoint'].iloc[peaks], data['qPCR'].iloc[peaks], marker='o', linestyle='', color='red', label='Peaks')
plt.xlabel('Timepoint (weeks)')
plt.ylabel('Parasitemia (Log-transformed qPCR)')
plt.legend()
plt.show()


## Third trial (Defining threshold using the mean and standard deviation of multiple random samplings of 
# the data)
# Define the number of random samples
num_samples = 100

# Define the sample size
sample_size = 50

# Create an empty list to store the means for each timepoint
all_means = []

# Perform random sampling and calculate means for each timepoint
for _ in range(num_samples):
    # Randomly choose sample_size number of kids
    random_sample = first_qpcr['Kid'].drop_duplicates().sample(n=sample_size, random_state=None)
    
    # Filter the data for the randomly chosen kids
    random_sample_data = first_qpcr[first_qpcr['Kid'].isin(random_sample)]
    
    # Calculate the mean qPCR value for each timepoint
    means = random_sample_data.groupby('Timepoint')['qPCR'].mean()
    
    # Append the means to the list
    all_means.append(means)

# Convert the list of means to a DataFrame
all_means_df = pd.DataFrame(all_means).T

# Plot all the resulting curves
plt.figure(figsize=(12, 6))
for i in range(num_samples):
    plt.plot(all_means_df.index, all_means_df.iloc[:, i], linestyle='-', alpha=0.3)

plt.xlabel('Timepoint (weeks)')
plt.ylabel('Mean Parasitemia (Log-transformed qPCR)')
plt.show()