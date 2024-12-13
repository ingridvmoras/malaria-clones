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
    def identify_peaks(self, df, col,lod=0, id):
        results = []
        for individual, group in df.groupby(id):
            if group.empty or len(group) == 1:
                continue
            
            #Using findpeaks library to identify peaks 
            #edogran/findpeaks: A Python library for finding peaks in 1D data. (github.com)
            #https://github.com/erdogant/findpeaks/
            
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
        std_data = df.groupby(id)[col].std()
        min, max = std_data.quantile([0.25, 0.75])
        
        for kid, group in df.groupby(id):
            height = lod
            prominence = 0
            threshold = 0.06
            
            std = group[col].std()
            if std < min:
                height = lod
                threshold -= 0.01
            elif std > max:
                height = lod + 1
                threshold += 0.01
                prominence += 1
                
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
                    'std': std
                }
                
                if (row_data['log2FoldChange'] < 1 and row_data['log2FoldChange'] == result['Prominence'] and 
                    result['RightBase'] - result['LeftBase'] == 2):
                    result['peak'] = False
                    result['falsetype'] = 'type1'
                elif (row_data['log2FoldChange'] < 1 and result['Prominence'] == result['RightThreshold']):
                    result['peak'] = False
                    result['falsetype'] = 'type2'
                elif (row_data['log2FoldChange'] < 1):
                    result['peak'] = False
                    result['falsetype'] = 'type3'
                else:
                    result['peak'] = True

                results.append(result)
        
        results = pd.DataFrame(results) 
        return results

class peak_finding:
    @staticmethod
    def method(method: str) -> PeakIdentifier:
        if method == 'topology':
            return TopologyPeakIdentifier()
        elif method == 'local':
            return LocalMaximaPeakIdentifier()
        else:
            raise ValueError(f"Unknown method: {method}")


def identify_peaks(df, col, lod, id, method):
    peak_identifier = peak_finding.method(method)
    return peak_identifier.identify_peaks(df, col, lod, id)