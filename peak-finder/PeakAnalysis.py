import scipy.stats as stats

class PeakAnalysis:
    def check_normality(self, data):
        mean_fold_change = data.groupby('Kid')['log2_qPCR'].std()
        k2, p = stats.normaltest(mean_fold_change)
        if p < 0.05:
            return 'Not Normal'
        else:
            return 'Normal'
