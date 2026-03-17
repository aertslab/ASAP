import os
import pycisTopic
import pandas as pd
import pickle
import gzip
import argparse
from pycisTopic.lda_models import *

print('importing done')
print('defining functions and classes....')

## Arguments
def make_argument_parser():
    """
    Creates an ArgumentParser to read the options for this script from
    sys.argv
    """
    parser = argparse.ArgumentParser(
        description="Load assigned topics and save topic model",)

    parser.add_argument('--cistopic_obj', '-cto', type=str, required=True,
                        help='cistopic object')

    parser.add_argument('--mallet', '-m', type=str, required=True,
                        help='Path to mallet.')

    parser.add_argument('--output_file', '-out', type=str, required=True,
                        help='output_file')

    return parser

def main():
    """
    The main executable function
    """
    print('make argument parser')
    parser = make_argument_parser()
    args = parser.parse_args()
    
    ## Print the arguments
    cto_path = args.cistopic_obj
    print('cistopic_obj: ', cto_path)
    
    mallet_path = args.mallet
    print('mallet path: ', mallet_path)

    outfile = args.output_file
    print('output: ', outfile)

    os.environ['MALLET_MEMORY'] = '1750G'
    
    print('Opening cto....')
    infile = open(cto_path, 'rb')
    cistopic_obj = pickle.load(infile)
    infile.close()
    print('Opening cto done')
    print('creating corpus.....')
    LDAMallet.convert_binary_matrix_to_mallet_corpus_file(cistopic_obj.binary_matrix, outfile, mallet_path)
    print('creating corpus finished')
    print('-------------------------')
    print('JOB DONE!')
    print('-------------------------')
  

if __name__ == "__main__":
    main()