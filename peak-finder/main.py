from PeakDetectionPipeline import PeakDetectionPipeline

if __name__ == "__main__":
    data_filepath = '..\\data\\dataset_Kalifabougou.csv'
    pipeline = PeakDetectionPipeline(data_filepath)
    pipeline.run()
