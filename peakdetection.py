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
        
        if 'Method' not in df.columns:
            df['Method'] = 'TPH'
            
        df['peak'] = False
        df['valley']= False 
        
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
            
            peak_values = peaks['df']['y'].values
            peak_indices = group[group[col].isin(peak_values)].index
            true_peak_indices = peak_indices[peaks['df']['peak'].values]
            #df['Method'] = df['Method'].where(~df.index.isin(true_peak_indices), 'TPH') #Topology Persistent Homology method
            df.loc[peak_indices, 'peak'] = peaks['df']['peak'].values
            df.loc[peak_indices, 'score'] = peaks['df']['score'].values
            df.loc[peak_indices, 'valley'] = peaks['df']['valley'].values
    
        return df
    
class LocalMaximaPeakIdentifier(PeakIdentifier):
    def identify_peaks(self, df, **kwargs):
        id = kwargs.get('id') # Individual identifier
        col = kwargs.get('col') # Column to find peaks in
        lod = kwargs.get('lod') # Level of detection
        win = kwargs.get('win') # Window size or distance
        
        std_data = df.groupby(id)[col].std()
        min, max = std_data.quantile([0.25, 0.75])
        
        df['peak'] = False
        df['height'] = np.nan
        df['Prominence'] = np.nan
        df['LeftBase'] = np.nan
        df['RightBase'] = np.nan
        df['LeftThreshold'] = np.nan
        df['RightThreshold'] = np.nan
        df['std'] = np.nan
        df['valley'] = False
        df['falsetype'] = np.nan
        df['Method'] = 'LM' #Local Maxima method
        # Autonomous parameter tuning based on standard deviation of the data 
        #Gupta, A., Gupta, S., Onumanyi, A.J. et al. A-TSPD: autonomous-two stage algorithm for robust peak 
        # detection in online time series. Cluster Comput 27, 4063–4076 (2024).
        # https://doi.org/10.1007/s10586-024-04369-8
        for individual, group in df.groupby(id):
            height = lod
            prominence = 0
            threshold = 0.06 
            
            std = group[col].std()
            if std < min:
                height = lod
                threshold =- 0.01
                
            elif std > max:
                height = lod + 1
                threshold += 0.01
                prominence += 1
                
        
            peaks, properties = find_peaks(group[col], height=height, prominence=prominence, wlen=win, threshold=threshold)
            peak_index = group.iloc[peaks].index
            
            for i, idx in enumerate(peak_index):
                df.loc[idx, 'peak'] = True
                df.loc[idx, 'height'] = height
                df.loc[idx, 'PeakHeight'] = properties['peak_heights'][i]
                df.loc[idx, 'Prominence'] = properties['prominences'][i]
                df.loc[idx, 'LeftBase'] = properties['left_bases'][i]
                df.loc[idx, 'RightBase'] = properties['right_bases'][i]
                df.loc[idx, 'LeftThreshold'] = properties['left_thresholds'][i]
                df.loc[idx, 'RightThreshold'] = properties['right_thresholds'][i]
                df.loc[idx, 'std'] = std
                #df.loc[idx, 'Method'] = 'LM' #Local Maxima method
                
                # Identification of 'false' peaks useful for quality control
                if (group.loc[idx, 'log2FoldChange'] < 1 and group.loc[idx, 'log2FoldChange'] == properties['prominences'][i] and 
                    properties['right_bases'][i] - properties['left_bases'][i] == 2):
                    df.loc[idx, 'peak'] = False
                    df.loc[idx, 'falsetype'] = 'type1'
                    df.loc[idx, 'valley'] = True
                elif (group.loc[idx, 'log2FoldChange'] < 1 and properties['prominences'][i] == properties['right_thresholds'][i]):
                    df.loc[idx, 'peak'] = False
                    df.loc[idx, 'falsetype'] = 'type2'
                    df.loc[idx, 'valley'] = True
                elif (group.loc[idx, 'log2FoldChange'] < 1):
                    df.loc[idx, 'peak'] = False
                    df.loc[idx, 'falsetype'] = 'type3'
                    df.loc[idx, 'valley'] = True
                else:
                    df.loc[idx, 'peak'] = True
        
        return df


class S1PeakIdentifier(PeakIdentifier):
    
    #Implementation of S1 peak finding algorithm described by 
    # Palshikar, Girish. (2009). Simple Algorithms for Peak Detection in Time-Series. 
    def identify_peaks(self, df, **kwargs):
        id = kwargs.get('id')
        col = kwargs.get('col')
        win = kwargs.get('win')
        
        df['peak'] = False
        df['Method'] = 'S1' #S1 method
        for individual, group in df.groupby(id):
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
             #df.loc[group.index[max_list], 'Method'] = 'S1'
        return df  




class peak_finding:
    @staticmethod
    def method(method: str) -> PeakIdentifier:
        if method == 'TPH':
            return TopologyPeakIdentifier()
        elif method == 'LM':
            return LocalMaximaPeakIdentifier()
        elif method == 'S1':
            return S1PeakIdentifier()
        else:
            raise ValueError(f"Unknown method: {method}")
        
        
#Function to identify peaks

def identify_peaks(df, method, **kwargs):
    """_summary_

    Args:
        df (_df_): Dataset with parasitemia time-series data. 
        method (_str_): Method to use for peak identification ("TPH', 'LM', 'S1')

    Returns:
        _type_: DataFrame with identified peaks
    """    
    peak_identifier = peak_finding.method(method)
    return peak_identifier.identify_peaks(df, **kwargs)