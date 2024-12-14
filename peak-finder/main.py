from PeakDetectionPipeline import PeakDetectionPipeline
from utils import dataset_statistics, joinSequences, createFileSequences, embeddingsGeneration
from dataloader import load_data

if __name__ == "__main__":
    data_filepath = '..\\data\\dataset_Kalifabougou.csv'
    pipeline = PeakDetectionPipeline(data_filepath)
    pipeline.run()
