from pipeline import PeakDetectionPipeline


if __name__ == "__main__":
    data_filepath = 'C:\\Users\\isabe\\OneDrive\\Escritorio\\malaria-clones\\peak-finder\\data\\dataset_Kalifabougou.csv' #Replace with the path to your data file
    pipeline = PeakDetectionPipeline(data_filepath)
    pipeline.run()
 