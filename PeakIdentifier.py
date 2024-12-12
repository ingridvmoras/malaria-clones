from __future__ import annotations
from abc import ABC, abstractmethod
from scipy.signal import find_peaks
from findpeaks import findpeaks
import pandas as pd
import numpy as np

class PeakIdentifier(ABC):
    @abstractmethod
    def identify_peaks(self, df, col, lod, id):
        pass

class TopologyPeakIdentifier(PeakIdentifier):
    def identify_peaks(self, df, col, lod, id):
        results = []
        for individual, group in df.groupby(id):
            if group.empty or len(group) == 1:
                continue
            fp = findpeaks(method='topology', lookahead=1)
            peaks = fp.peaks1d(X=group[col], method='topology')
            if peaks is None or 'df' not in peaks or peaks['df'] is None or len(peaks['df']) == 0:
                continue
            peaks['df'].rename(columns={'y': col}, inplace=True)
            peaks['df']['Timepoint'] = group['Timepoint'].values
            peaks['df'][id] = individual
            peaks['df']['peak'] = peaks['df']['peak'].astype(bool)
            peaks['df']['valley'] = peaks['df']['valley'].astype(bool)
            results.append(peaks['df'])
        final_df = pd.concat(results, ignore_index=True)
        return final_df

class LocalMaximaPeakIdentifier(PeakIdentifier):
    def identify_peaks(self, df, col, lod, id):
        results = []
        for kid, group in df.groupby(id):
            height = lod
            prominence = 0
            threshold = 0.06
            std = group[col].std()
            peaks, properties = find_peaks(group[col], height=height, prominence=prominence, wlen=5, threshold=threshold)
            peak_rows = group.iloc[peaks]
            for i, (idx, row_data) in enumerate(peak_rows.iterrows()):
                result = {
                    'Kid': kid,
                    'height': height,
                    'Timepoint': row_data['Timepoint'],
                    col: row_data[col],
                    'PeakHeight': properties['peak_heights'][i],
                    'Prominence': properties['prominences'][i],
                    'log2FoldChange': row_data['log2FoldChange'],
                    'LeftBase': properties['left_bases'][i],
                    'RightBase': properties['right_bases'][i],
                    'LeftThreshold': properties['left_thresholds'][i],
                    'RightThreshold': properties['right_thresholds'][i],
                    'std': std,
                    'peak': True
                }
                results.append(result)
        results = pd.DataFrame(results)
        return results