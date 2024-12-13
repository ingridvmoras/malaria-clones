import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

# set the directory containing the CSV files
directory = '../outcome/'
all_files = glob.glob(os.path.join(directory, "infection_duration_analysis*.csv"))
current_date = datetime.now().strftime("%Y-%m-%d")

for file in all_files:
    max_zeroes = file.split('_')[-3]
    data = pd.read_csv(file)
    count = data.groupby(['Timepoint', 'first_appearance']).size().reset_index(name='count')
    count = count[count['first_appearance'] == 1]
    count['Timepoint'] = count['Timepoint'].replace(44, "May13")
    
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Timepoint', y='count', data=count, color='#408F8A')
    plt.xlabel('Timepoint')
    plt.ylabel('First appearances per timepoint')
    plt.xticks(rotation=0, fontsize=10)
    plt.yticks(fontsize=10)
    plt.title(f'First Appearances Per Timepoint (max_zeroes={max_zeroes})')
    plot_filename = f"first_appearances_per_timepoint_{max_zeroes}_zeroes_{current_date}.png"
    plt.savefig(os.path.join('../plots/', plot_filename))
    plt.close()
    
    myTab = pd.crosstab(data['first_appearance'], data['peak'])
    real1 = data[data['first_TP'] == True]
    myTab2 = pd.crosstab(real1['first_appearance'], real1['peak'])


    FA_real1 = real1.groupby('first_appearance').size().reset_index(name='count')  
    FA_real = data.groupby('first_appearance').size().reset_index(name='count')  

    first_appearance_per_100_peaks_observed = round(myTab2.loc[1, 1] / (myTab2.loc[1, 1] + myTab2.loc[1, 0]) * 100) #17
    
    

