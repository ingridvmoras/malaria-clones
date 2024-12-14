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
        df['method'] = 'topology'
        results = []
        for individual, group in df.groupby(id):
        
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
                df.loc[peak_indices, 'peak'] = True
                df.loc[peak_indices, 'score'] = peaks['df']['score'].values
                df.loc[peak_indices, 'valley'] = peaks['df']['valley'].values
    
        return df
    
class LocalMaximaPeakIdentifier(PeakIdentifier):
    def identify_peaks(self, df, **kwargs):
        id = kwargs.get('id') # Individual identifier
        col = kwargs.get('col') # Column to find peaks in
        lod = kwargs.get('lod') # Level of detection
        win = kwargs.get('win') # Window size or distance
        results = []
        std_data = df.groupby(id)[col].std()
        min, max = std_data.quantile([0.25, 0.75])
        
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
                
        
            peaks, properties = find_peaks(group)[col], height=height, prominence=prominence, wlen=win, threshold=threshold)
            peak_rows = group.iloc[peaks]
            
            for i, (idx, row_data) in enumerate(peak_rows.iterrows()):
                
                result = {
                    f'{id}': individual,
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
                    'valley': False
                }
                
                
                #Identification of 'false' peaks useful for quality control
                if (row_data['log2FoldChange'] < 1 and row_data['log2FoldChange'] == result['Prominence'] and 
                    result['RightBase'] - result['LeftBase'] == 2):
                    result['peak'] = False
                    result['falsetype'] = 'type1'
                    result['valley'] = True
                elif (row_data['log2FoldChange'] < 1 and result['Prominence'] == result['RightThreshold']):
                    result['peak'] = False
                    result['falsetype'] = 'type2'
                    result['valley'] = True
                elif (row_data['log2FoldChange'] < 1):
                    result['peak'] = False
                    result['falsetype'] = 'type3'
                    result['valley'] = True
                else:
                    result['peak'] = True

                results.append(result)
        
        results = pd.DataFrame(results) 
        return results


class S1PeakIdentifier(PeakIdentifier):
    
    #Implementation of S1 peak finding algorithm described by 
    # Palshikar, Girish. (2009). Simple Algorithms for Peak Detection in Time-Series. 
    def identify_peaks(self, df, **kwargs):
        id = kwargs.get('id')
        col = kwargs.get('col')
        win = kwargs.get('win')
        
        df['peak'] = False

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
    """_summary_

    Args:
        df (_type_): Dataset with parasitemia time-series data. 
        method (_type_): Method to use for peak identification ("topology', 'local', 's1')

    Returns:
        _type_: DataFrame with identified peaks
    """    
    peak_identifier = peak_finding.method(method)
    return peak_identifier.identify_peaks(df, **kwargs)