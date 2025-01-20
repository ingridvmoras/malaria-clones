import matplotlib.pyplot as plt
import seaborn as sns
# import seaborn.objects as so
import pandas as pd
import utils 
import os
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import ConfusionMatrixDisplay
import numpy as np


sns.set_style("white") 
cb_palette = ["#999999", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
    # http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
    # http://jfly.iam.u-tokyo.ac.jp/color/


def plot_mean_std(data, col, id):
    mean_std = data.groupby(id)[col].std()
    plt.figure(figsize=(10, 6))
    plt.hist(mean_std, bins=30, edgecolor='black')
    plt.axvline(mean_std.mean(), color='red', linestyle='dotted', linewidth=2, label=f'Mean: {mean_std.mean():.2f}')
    quantiles = mean_std.quantile([0.25, 0.5, 0.75])
    for q in quantiles.index:
        plt.axvline(quantiles[q], color='blue', linestyle='dotted', linewidth=2, label=f'{int(q*100)}th percentile: {quantiles[q]:.2f}')
    plt.xlabel('Mean Standard Deviation (log2)')
    plt.ylabel('Frequency')
    plt.show()
    
def plot_methods(df, id, col, f_name='methods_plot', save_path='plots', num_sample=1):
    os.makedirs(save_path, exist_ok=True)
    sample = utils.sample_random(df, id, num_sample)
    timepoints = sorted(df['Timepoint'].unique())
    pdf_path = os.path.join(save_path, f'{f_name}.pdf')
    with PdfPages(pdf_path) as pdf:
        for unique_id in sample[id].unique():
            sample_id = sample[sample[id] == unique_id]
            plt.figure(figsize=(18, 8), dpi=300)
            
            methods = sample_id['Method'].dropna().unique()
            num_methods = len(methods)
            method_colors = {method: cb_palette[i % len(cb_palette)] for i, method in enumerate(methods)}
            
            # Fig A
            for i, method in enumerate(methods):
                ax = plt.subplot(num_methods, 2, 2 * i + 1)
                method_data = sample_id[sample_id['Method'] == method]
                plot_data = sample_id.drop_duplicates(subset=['Timepoint', 'Kid'], keep='first').copy()
                sns.lineplot(data=plot_data, x='Timepoint', y=col, ax=ax, label=method, color=method_colors[method])
                sns.scatterplot(data=method_data[method_data['peak']], x='Timepoint', y=col, ax=ax, color='red', label='Peaks')
                ax.set_xticks(timepoints)
                ax.set_xlabel('')
                ax.set_ylabel(f'{col}(p/uL)')
                if i == num_methods - 1:
                    ax.set_xlabel('Timepoint')
                ax.legend(frameon=True, fancybox=True, framealpha=1, borderpad=0)
                ax.text(-0.12, 0.5, f'{method}'.upper(), transform=ax.transAxes, fontsize=12, fontweight='bold', va='center', ha='center', rotation=90)
            
            pdf.savefig(bbox_inches='tight')
            plt.close()
            
    # Fig B
    plt.figure(figsize=(8, 6), dpi=300)
    methods_order = ['TPH', 'LM', 'S1']
    peak_counts = df[df['peak'] == True].groupby(['Method', id]).size().reset_index(name='peak_count')
    peak_counts['Method'] = pd.Categorical(peak_counts['Method'], categories=methods_order, ordered=True)
    sns.swarmplot(data=peak_counts, x='Method', y='peak_count', palette=cb_palette, hue='Method', order=methods_order)
    plt.xlabel('Method')
    plt.ylabel(f'Number of Peaks per {id}')
    plt.tight_layout()
    plt.yticks(np.arange(0, peak_counts['peak_count'].max() + 1, step=1))
    plt.savefig(f'{save_path}/swarmplot_{f_name}', bbox_inches='tight')
    plt.close()
        

def plot_matrix(method, pdf, metrics, conf_matrix, contingency_table=None):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    
    cm_display = ConfusionMatrixDisplay(conf_matrix, display_labels = [False, True])
    cm_display.plot(ax=axes[0], cmap='Blues') 
    axes[0].set_xlabel('Predicted')
    axes[0].set_ylabel('Ground Truth') 
    axes[0].set_title(f'Confusion Matrix for {method}')
    
    if contingency_table is not None:
        sns.heatmap(contingency_table, annot=True, fmt='d', cmap='Blues', 
                    xticklabels=['False', 'True'], 
                    yticklabels=['False', 'True'], ax=axes[1])
        axes[1].set_xlabel('Valley')
        axes[1].set_ylabel('Peak')
        axes[1].set_title(f'Contingency Matrix for {method}')
        
    plt.show()
    pdf.savefig(fig)
    plt.close(fig)
    
    print(f'Metrics for {method}:')
    for metric_name, metric_value in metrics.items():
        print(f'{metric_name}: {metric_value}')

def plot_simulations(all_data, plot_dir, current_date):
    methods = all_data['method'].unique()
    for method in methods:
        method_data = all_data[all_data['method'] == method]
        plt.figure(figsize=(10, 6))
        sns.histplot(data=method_data, x='number_of_first_appearances_100_peaks', hue='type', element='step', stat='count', binwidth=1)
        plt.title(f'First Appearances per 100 Peaks Observed - Method: {method}')
        plt.xlabel('Number of First Appearances per 100 Peaks')
        plt.ylabel('Count')
        plt.legend( labels=['Weighted', 'Unweighted'])
        plot_filename = f'first_appearances_per_100_peaks_{method}_{current_date}.png'
        plt.savefig(os.path.join(plot_dir, plot_filename))
        plt.close()

def plot_first_appearances(df, plot_dir):
    sns.lineplot(x='max_zeroes', y='first_appearance_per_100_peaks_observed', hue='type', style='method', markers=True, dashes=False, data=df)
    plt.xlabel('Skips')
    plt.ylabel('First Appearances per 100 Peaks Observed')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(ticks=sorted(df['max_zeroes'].unique()))  
    plt.tight_layout()
    plot_filename = f'first_appearances_per_100_peaks.png'
    plt.savefig(os.path.join(plot_dir, plot_filename))
    plt.close()
