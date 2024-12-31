import scipy.stats as stats
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages 
import plots as p

class PeakAnalysis:
    def check_normality(self, data, id, col):
        mean_fold_change = data.groupby(id)[col].std()
        k2, p = stats.normaltest(mean_fold_change)
        print(f"Test statistic: {k2}")  # Accessing k2
        if p < 0.05:
            return 'Not Normal'
        else:
            return 'Normal'
    def analyze_methods(self, df, id, path='outcome\\method_analysis.pdf'):

        df = df.dropna(subset=['Method'])
        df.set_index(['Method','Timepoint',id], append=True, inplace=True)
        df['peak_count'] = df.groupby(level=[id, 'Timepoint'])['peak'].transform('sum')
        df = df.reset_index()
        
        df['ground_truth'] = (df['peak_count'] == len(df['Method'].dropna().unique())).astype(int)
        unique_methods = df['Method'].unique()
        
        with PdfPages(path) as pdf:
            for method in unique_methods:
                df_method = df[df['Method'] == method]
                y_true = df_method['ground_truth'].astype(int)
                y_pred = df_method['peak'].astype(int)
                
                metrics = {
                    'Accuracy': accuracy_score(y_true, y_pred),
                    'Precision': precision_score(y_true, y_pred),
                    'Recall': recall_score(y_true, y_pred),
                    'F1 Score': f1_score(y_true, y_pred)
                }
                conf_matrix = confusion_matrix(y_true, y_pred)
                
                contingency_table = None
                if method in ['TPH', 'LM']:
                    contingency_table = pd.crosstab(df_method['peak'], df_method['valley'])
                
                p.plot_matrix(method, pdf, metrics, conf_matrix, contingency_table)

