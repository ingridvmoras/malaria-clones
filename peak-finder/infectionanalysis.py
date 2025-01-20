import pandas as pd
import glob
import os
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

#Code modified from apperances_Sept2023.R by Manuela Carrasquilla

def analyze_infection(directory='outcome',files="infection_duration_analysis*.csv",plot_dir='plots'):
    all_files = glob.glob(os.path.join(directory,files))
    current_date = datetime.now().strftime("%Y-%m-%d")

    for file in all_files:
        info = []
        max_zeroes = file.split('_')[-4]
        method = file.split('_')[-1].split('.')[0]
        data = pd.read_csv(file)
        count = data.groupby(['Timepoint', 'first_appearance']).size().reset_index(name='count')
        count = count[count['first_appearance'] == 1]
        count['Timepoint'] = count['Timepoint'].replace(13, "May13")

        plt.figure(figsize=(10, 6))
        sns.barplot(x='Timepoint', y='count', data=count, color='#408F8A')
        plt.xlabel('Timepoint')
        plt.ylabel('First appearances per timepoint')
        plt.xticks(rotation=0, fontsize=10)
        plt.yticks(fontsize=10)
        plt.title(f'First Appearances Per Timepoint (max_zeroes={max_zeroes})')
        plot_filename = f"first_appearances_per_timepoint_{max_zeroes}_zeroes_{current_date}.png"
        plt.savefig(os.path.join(plot_dir, plot_filename))
        plt.close()

        myTab = pd.crosstab(data['first_appearance'], data['peak'])
        print(myTab)
        sns.heatmap(myTab, annot=True, fmt='d', cmap='Blues', 
                    xticklabels=['False', 'True'], 
                    yticklabels=['False', 'True'])
        plt.xlabel('Peak')
        plt.ylabel('First appearance')
        plot_filename = f"crosstab_{max_zeroes}_zeroes_{method}_{current_date}.png"
        plt.savefig(os.path.join(plot_dir, plot_filename))
        plt.close()

        FA_real = data.groupby('first_appearance').size().reset_index(name='count')
        print(f'Number of first appearances: {FA_real.loc[1, "count"]}')

        first_appearance_per_100_peaks_observed = round(myTab.loc[1, 1] / (myTab.loc[0, 1] + myTab.loc[1, 1]) * 100)
        print(f'Number of first appearances per 100 peaks observed: {first_appearance_per_100_peaks_observed}')

        types = data['Type'].unique()
        for t in types:
            subset = data[data['Type'] == t]
            myTab3 = pd.crosstab(subset['first_appearance'], subset['peak'])
            first_appearance_per_100_peaks_observed = round(myTab.loc[1, 1] / (myTab.loc[0, 1] + myTab.loc[1, 1]) * 100)
            info.append([method, max_zeroes, t, first_appearance_per_100_peaks_observed])

        info_df = pd.DataFrame(info, columns=['method', 'max_zeroes', 'type', 'first_appearance_per_100_peaks_observed'])
        info_df.to_csv(f'{directory}first_appearances_per_100_peaks_observed_{current_date}.csv', index=False)

def analyze_simulations(directory='plots', files='simulation_*'):
    all_files = glob.glob(os.path.join(directory, files))
    current_date = datetime.now().strftime("%Y-%m-%d")

    results = []

    for file in all_files:
        data = pd.read_csv(file)
        method = file.split('_')[-1].split('.')[0]
        data['method'] = method

        if 'weighted' in file:
            simulations = [group for _, group in data.groupby('simulation')]
            first_appearance_per_100_peaks = []
            simulation_number = []

            for k, simulation in enumerate(simulations, start=1):
                mytab = pd.crosstab(simulation['first_appearance'], simulation['is_peak'])
                simulation_number.append(k)
                first_appearance_per_100_peaks.append(round(mytab.loc[1, 1] / (mytab.loc[0, 1] + mytab.loc[1, 1]) * 100))

            weighted = pd.DataFrame({'simulation': simulation_number, 'number_of_first_appearances_100_peaks': first_appearance_per_100_peaks})
            weighted['type'] = 'weighted'
            weighted['method'] = method
            results.append(weighted)
        
        elif 'unweighted' in file:
            simulations = [group for _, group in data.groupby('simulation')]
            first_appearance_per_100_peaks = []
            simulation_number = []

            for k, simulation in enumerate(simulations, start=1):
                mytab = pd.crosstab(simulation['first_appearance'], simulation['is_peak'])
                simulation_number.append(k)
                first_appearance_per_100_peaks.append(round(mytab.loc[1, 1] / (mytab.loc[0, 1] + mytab.loc[1, 1]) * 100))

            unweighted = pd.DataFrame({'simulation': simulation_number, 'number_of_first_appearances_100_peaks': first_appearance_per_100_peaks})
            unweighted['type'] = 'unweighted'
            unweighted['method'] = method
            results.append(unweighted)
        
    all_data = pd.concat(results, ignore_index=True)
    methods = all_data['method'].unique()
    for method in methods:
        method_data = all_data[all_data['method'] == method]
        plt.figure(figsize=(10, 6))
        sns.histplot(data=method_data, x='number_of_first_appearances_100_peaks', hue='type', element='step', stat='count', binwidth=1, kde=True)
        plt.title(f'First Appearances per 100 Peaks Observed - Method: {method}')
        plt.xlabel('Number of First Appearances per 100 Peaks')
        plt.ylabel('Count')
        plt.legend(title='Type')
        plt.savefig(f'{directory}/first_appearances_per_100_peaks_simulations_{method}_{current_date}.png')
        plt.close()

    return all_data


