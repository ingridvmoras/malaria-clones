from pipeline import PeakDetectionPipeline


if __name__ == "__main__":
    data_filepath = 'N:\\Mora\\project\\malaria-clones\\peak-finder\\data\\dataset_Kalifabougou.csv'
    pipeline = PeakDetectionPipeline(data_filepath)
    pipeline.run()
