from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from peak_finder.dataprocessor import DataProcessor
from peak_finder.peakdetection import identify_peaks
from peak_finder.peakanalysis import PeakAnalysis
from peak_finder.infectionanalysis import analyze_infection, analyze_simulations
import peak_finder.plots as plots
import peak_finder.utils as f
import os


class PeakDetectionPipeline:
    def __init__(self, data_filepath, **kwargs):
        self.data_processor = DataProcessor(data_filepath)
        self.kwargs = kwargs
        self.peak_analysis = PeakAnalysis()
        self.data_filepath = data_filepath
        

    def run(self):
        '''
        Executes the main pipeline for peak detection and analysis.
        This method performs the following steps:
        1. Creates necessary directories for output and plots.
        2. Loads and preprocesses the data.
        3. Identifies peaks using different methods (TPH, LM, S1).
        4. Merges results from different peak detection methods.
        5. Analyzes the performance of peak detection methods.
        6. Plots the results.
        7. Saves the results to CSV files with peak data merged with original data.
        8. Runs simulations and saves the results.
        9. Analyzes the number of first appearances per 100 peaks observed for each method.
        Raises:
            Exception: If any step in the pipeline fails.
        '''    
        output_dir = 'outcome'
        plots_dir = 'plots'
        os.makedirs(plots_dir, exist_ok=True)
        os.makedirs(output_dir, exist_ok=True)
        
        # Load and preprocess data 
        data = self.data_processor.load_data(self.data_filepath) # Load data
        
        if self.kwargs.get('preprocess', True):
            preprocessed_data, original_data = self.data_processor.preprocess_data(data) # Preprocess data if needed
            processed_data, num_ids= self.data_processor.process_data(preprocessed_data)
        else:
            processed_data, num_ids = self.data_processor.process_data(data)
           
        # Identify peaks using different methods
        
        
        
        df_tph = identify_peaks(processed_data.copy(), method='TPH', **self.kwargs) #Topology Persistent Homology
        df_lm = identify_peaks(processed_data.copy(), method='LM', **self.kwargs) #Local Maxima
        df_s1 = identify_peaks(processed_data.copy(), method='S1', **self.kwargs) #S1
        
        # Merge results from different peak detection methods
        mdf = f.mergedf([df_tph, df_lm, df_s1], col=self.kwargs['col'], id=self.kwargs['id'])

        
        mdf.to_csv(f'{output_dir}/peaks.csv', index=False) #Returns only peak data 
    
        # Analyze performance of peak detection methods
        self.peak_analysis.analyze_methods(mdf, id= self.kwargs['id'],output_dir=output_dir)
        
        # Plot results 
        plots.venndiagram(df_tph,df_lm,df_s1,self.kwargs['id'],save_path=plots_dir)
        plots.plot_methods(mdf, id=self.kwargs['id'], col=self.kwargs['col'], f_name='methods', save_path=plots_dir, num_sample=num_ids)
        
        # Save results to csv with peak data merged with original data
        methods = {'TPH': df_tph, 'LM': df_lm, 'S1': df_s1}
        for method, df_peak in methods.items():
            df_peak['peak'] = df_peak['peak'].astype(int)
            df = f.merge_peak_data(original_data, df_peak, method=method)
            f.first_appearances(df, id=self.kwargs['id'], method=method, output_dir=output_dir)
            s1= f.run_simulation(df, id=self.kwargs['id'], rounds=1000, weighted=True)
            s2= f.run_simulation(df, id=self.kwargs['id'], rounds=1000, weighted=False)
            s1.to_csv(f'{output_dir}/simulation_weighted_{method}.csv', index=False)
            s2.to_csv(f'{output_dir}/simulation_unweighted_{method}.csv', index=False)
        
        #Get the number of first apparitions per 100 peaks observed for each method 
        
        analyze_infection(output_dir,files="infection_duration_analysis*.csv", plot_dir=plots_dir)
        data= analyze_simulations(output_dir,files='simulation_*', plot_dir=plots_dir)
        

