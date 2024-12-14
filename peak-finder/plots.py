import matplotlib.pyplot as plt
import seaborn as sns
import seaborn.objects as so


sns.set_style("ticks") 
cb_palette = ["#999999", "#E69F00", "#56B4E9", "#009E73",
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
    # http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
    # http://jfly.iam.u-tokyo.ac.jp/color/

class plots:
    def plot_mean_fold_change(self, data):
        mean_fold_change = data.groupby('Kid')['log2_qPCR'].std()
        plt.figure(figsize=(10, 6))
        plt.hist(mean_fold_change, bins=30, edgecolor='black')
        plt.axvline(mean_fold_change.mean(), color='red', linestyle='dotted', linewidth=2, label=f'Mean: {mean_fold_change.mean():.2f}')
        quantiles = mean_fold_change.quantile([0.25, 0.5, 0.75])
        for q in quantiles.index:
            plt.axvline(quantiles[q], color='blue', linestyle='dotted', linewidth=2, label=f'{int(q*100)}th percentile: {quantiles[q]:.2f}')
        plt.xlabel('Mean Fold Change (log2)')
        plt.ylabel('Frequency')
        plt.show()
