from __future__ import annotations
from abc import ABC, abstractmethod
from scipy.signal import find_peaks
from findpeaks import findpeaks
import pandas as pd
import numpy as np

class PeakIdentifier(ABC):
    @abstractmethod
    def identify_peaks(self, df, **kwargs):
        pass

class TopologyPeakIdentifier(PeakIdentifier):
    def identify_peaks(self, df, **kwargs):
        id = kwargs.get('id')
        col = kwargs.get('col')
        
        results = []
        for individual, group in df.groupby(id):
            if group.empty or len(group) == 1:
                continue
            
            # Using findpeaks library to identify peaks 
            # Taskesen, E. (2020). findpeaks is for the detection of peaks
            # and valleys in a 1D vector and 2D array (image). (Version 2.3.1) 
            # [Computer software]. https://erdogant.github.io/findpeaks
            
            fp = findpeaks(method='topology', lookahead=1)
            peaks = fp.peaks1d(X=group[col], method='topology')
            if peaks is None or 'df' not in peaks or peaks['df'] is None or len(peaks['df']) == 0:
                continue
            peaks['df'].rename(columns={'y': col}, inplace=True)
            peaks['df']['Timepoint'] = group['Timepoint'].values
            peaks['df'][id] = individual
            peaks['df']['peak'] = peaks['df']['peak'].astype(bool)
            results.append(peaks['df'])
        
        return pd.concat(results, ignore_index=True) if results else pd.DataFrame()

class LocalMaximaPeakIdentifier(PeakIdentifier):
    def identify_peaks(self, df, **kwargs):
        id = kwargs.get('id') # Individual identifier
        col = kwargs.get('col') # Column to find peaks in
        lod = kwargs.get('lod') # Level of detection
        win = kwargs.get('win') # Window size or distance
        
        results = []
        for individual, group in df.groupby(id):
            if group.empty or len(group) == 1:
                continue
            
            peaks, _ = find_peaks(group[col], height=lod, distance=win)
            if len(peaks) == 0:
                continue
            
            peak_df = group.iloc[peaks].copy()
            peak_df['peak'] = True
            results.append(peak_df)
        
        return pd.concat(results, ignore_index=True) if results else pd.DataFrame()

class S1PeakIdentifier(PeakIdentifier):
    
    #Implementation of S1 peak finding algorithm described by 
    # Palshikar, Girish. (2009). Simple Algorithms for Peak Detection in Time-Series. 
    def identify_peaks(self, df, **kwargs):
        id = kwargs.get('id')
        col = kwargs.get('col')
        win = kwargs.get('win')
        
        df['peak'] = False

        for kid, group in df.groupby(id):
             data = group[col].values
             data_extended = np.concatenate([np.zeros(win), data, np.zeros(win)])
             max_list = []

             for i, value in enumerate(data_extended):
                 if (i >= win) and (i < len(data_extended) - win):
                    try:
                        max_left = data_extended[(i - win):i + 1].max()
                        max_right = data_extended[i:(i + win) + 1].max()
                        chek_value = data_extended[i] - ((max_left + max_right) / 2)
                    except ValueError:
                        pass

                    if (chek_value >= 0):
                        max_list.append(i - win)

             df.loc[group.index[max_list], 'peak'] = True
        return df




class peak_finding:
    @staticmethod
    def method(method: str) -> PeakIdentifier:
        if method == 'topology':
            return TopologyPeakIdentifier()
        elif method == 'local':
            return LocalMaximaPeakIdentifier()
        elif method == 's1':
            return S1PeakIdentifier()
        else:
            raise ValueError(f"Unknown method: {method}")
        
        
#Function to identify peaks

def identify_peaks(df, method, **kwargs):
    peak_identifier = peak_finding.method(method)
    return peak_identifier.identify_peaks(df, **kwargs)