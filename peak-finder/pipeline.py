import numpy as np
import pandas as pd
import seaborn as sns
from dataprocessor import DataProcessor
from peakdetection import identify_peaks
from peakanalysis import PeakAnalysis
import plots
import utils as f
import os

class PeakDetectionPipeline:
    def __init__(self, data_filepath):
        self.data_processor = DataProcessor(data_filepath)
        self.peak_analysis = PeakAnalysis()
        self.data_filepath = data_filepath

    def run(self):
        output_dir = 'outcome'
        os.makedirs(output_dir, exist_ok=True)
        
        # Load and preprocess data 
        data = self.data_processor.load_data(self.data_filepath) # Load data
        preprocessed_data, original_data = self.data_processor.preprocess_data(data) #Preprocess data if needed
        processed_data = self.data_processor.process_data(preprocessed_data) # Process data if needed
           
        # Identify peaks using different methods
        kwargs = {
            'id': 'Kid',          # Column name with ID of individuals
            'col': 'log2_qPCR',   # Col name with parasitemia data 
            'lod': np.log2(100),  # Level of detection for LocalMaximaPeakIdentifier
            'win': 5             # Window size 
        }
        
        df_tph = identify_peaks(processed_data.copy(), method='TPH', **kwargs) #Topology Persistent Homology
        df_lm = identify_peaks(processed_data.copy(), method='LM', **kwargs) #Local Maxima
        df_s1 = identify_peaks(processed_data.copy(), method='S1', **kwargs) #S1
        
        # Merged results from different peak detection methods
        mdf = f.mergedf([df_tph, df_lm, df_s1],col= kwargs['col'],id= kwargs['id'])
        
        # Analyze performance of peak detection methods
        self.peak_analysis.analyze_methods(mdf, id= kwargs['id'])
        
        # Plot results 
        plots.plot_methods(mdf,col= kwargs['col'],id= kwargs['id'])
        
        # Save results to csv with peak data merged with original data
        methods = {'TPH': df_tph, 'LM': df_lm, 'S1': df_s1}
        for method, df_peak in methods.items():
            df = f.merge_peak_data(original_data, df_peak, method=method)
            f.first_appearances(df, id=kwargs['id'], method=method, output_dir=output_dir)
           
