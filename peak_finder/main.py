import argparse
import numpy as np
from peak_finder.pipeline import PeakDetectionPipeline
def main():
    parser = argparse.ArgumentParser(description='Run the peak detection pipeline.')
    parser.add_argument('data_filepath', type=str, help='Path to the data file')
    parser.add_argument('--id', type=str, default='Kid', help='Column name with ID of individuals')
    parser.add_argument('--col', type=str, default='log2_qPCR', help='Column name with parasitemia data')
    parser.add_argument('--lod', type=float, default=np.log2(100), help='Level of detection for LocalMaximaPeakIdentifier')
    parser.add_argument('--win_s1', type=int, default=3, help='Window size for S1 method')
    parser.add_argument('--win_lm', type=int, default=6, help='Window size for LM method')
    parser.add_argument('--preprocess', action='store_true', help='Flag to indicate if preprocessing should be done')
    parser.add_argument('--falsetype', action='store_true', help='Flag to indicate if identification of false peak types should be done when using LM method.')
    
    args = parser.parse_args()

    kwargs = {
        'id': args.id,
        'col': args.col,
        'lod': args.lod,
        'win_s1': args.win_s1,
        'win_lm': args.win_lm,
        'preprocess': args.preprocess,
        'falsetype':args.falsetype
    }

    pipeline = PeakDetectionPipeline(args.data_filepath, **kwargs)
    pipeline.run()

if __name__ == "__main__":
    main()