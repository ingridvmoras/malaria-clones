import numpy as np  # Added this line
from DataProcessor import DataProcessor
from PeakIdentifier import TopologyPeakIdentifier, LocalMaximaPeakIdentifier
from PeakPlotter import PeakPlotter
from PeakAnalysis import PeakAnalysis

class PeakDetectionPipeline:
    def __init__(self, data_filepath):
        self.data_processor = DataProcessor()
        self.topology_identifier = TopologyPeakIdentifier()
        self.local_maxima_identifier = LocalMaximaPeakIdentifier()
        self.peak_plotter = PeakPlotter()
        self.peak_analysis = PeakAnalysis()
        self.data_filepath = data_filepath

    def run(self):
        data = self.data_processor.load_data(self.data_filepath)
        preprocessed_data = self.data_processor.preprocess_data(data)
        topology_peaks = self.topology_identifier.identify_peaks(preprocessed_data, 'log2_qPCR', np.log2(100), 'Kid')
        local_maxima_peaks = self.local_maxima_identifier.identify_peaks(preprocessed_data, 'log2_qPCR', np.log2(100), 'Kid')
        self.peak_plotter.plot_mean_fold_change(preprocessed_data)
        normality = self.peak_analysis.check_normality(preprocessed_data)
        print(f'Mean Fold Change Distribution: {normality}')
