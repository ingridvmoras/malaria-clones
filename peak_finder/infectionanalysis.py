import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from peak_finder.plots import plot_simulations, plot_first_appearances
from datetime import datetime

def analyze_infection(directory='outcome', files="infection_duration_analysis*.csv", plot_dir='plots'):
    all_files = glob.glob(os.path.join(directory, files))
    current_date = datetime.now().strftime("%Y-%m-%d")

    info = []  # Inicializar la lista fuera del bucle para acumular datos de todos los archivos

    for file in all_files:
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
        os.makedirs(plot_dir, exist_ok=True) 
        plt.savefig(os.path.join(plot_dir, plot_filename))
        plt.close()

        myTab = pd.crosstab(data['first_appearance'], data['peak'])
        print(f'{method} with max_zeroes={max_zeroes}')
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
            first_appearance_per_100_peaks_observed = round(myTab3.loc[1, 1] / (myTab3.loc[0, 1] + myTab3.loc[1, 1]) * 100)
            info.append([method, max_zeroes, t, first_appearance_per_100_peaks_observed])

    info_df = pd.DataFrame(info, columns=['method', 'max_zeroes', 'type', 'first_appearance_per_100_peaks_observed'])
    os.makedirs(directory, exist_ok=True)  
    info_df.to_csv(f'{directory}/first_appearances_per_100_peaks_observed_{current_date}.csv', index=False)
    plot_first_appearances(info_df, plot_dir)




def analyze_simulations(directory, files='simulation_*',plot_dir='plots'):
    all_files = glob.glob(os.path.join(directory, files))
    current_date = datetime.now().strftime("%Y-%m-%d")

    if not all_files:
        print(f"No files found in directory '{directory}' with pattern '{files}'")
        return

    results = []

    for file in all_files:
        print(f"Processing file: {file}")
        data = pd.read_csv(file)
        method = file.split('_')[-1].split('.')[0]
        data['method'] = method

        if 'weighted' in file:
            simulations = [group for _, group in data.groupby('simulation')]
            first_appearance_per_100_peaks = []
            simulation_number = []

            for k, simulation in enumerate(simulations, start=1):
                mytab = pd.crosstab(simulation['first_appearance'], simulation['peak'])
                simulation_number.append(k)
                first_appearance_per_100_peaks.append(round(mytab.loc[1, 1] / (mytab.loc[0, 1] + mytab.loc[1, 1]) * 100))

            weighted = pd.DataFrame({'simulation': simulation_number, 'number_of_first_appearances_100_peaks': first_appearance_per_100_peaks})
            weighted['type'] = 'weighted'
            weighted['method'] = method
            results.append(weighted)

        if 'unweighted' in file:
            simulations = [group for _, group in data.groupby('simulation')]
            first_appearance_per_100_peaks = []
            simulation_number = []

            for k, simulation in enumerate(simulations, start=1):
                mytab = pd.crosstab(simulation['first_appearance'], simulation['peak'])
                simulation_number.append(k)
                first_appearance_per_100_peaks.append(round(mytab.loc[1, 1] / (mytab.loc[0, 1] + mytab.loc[1, 1]) * 100))

            unweighted = pd.DataFrame({'simulation': simulation_number, 'number_of_first_appearances_100_peaks': first_appearance_per_100_peaks})
            unweighted['type'] = 'unweighted'
            unweighted['method'] = method
            results.append(unweighted)

    if not results:
        print("No results to concatenate")
        return

    all_data = pd.concat(results, ignore_index=True)

    plot_simulations(all_data, plot_dir, current_date)
    

    return all_data


