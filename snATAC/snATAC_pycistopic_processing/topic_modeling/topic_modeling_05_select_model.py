## Topic modeling evaluation: 
print('importing modules...')
import pandas as pd
import argparse
import pickle
import os
import re
import glob
import sys
import matplotlib.pyplot as plt
import plotnine as pn
from plotnine import *
from pycisTopic.cistopic_class import *
from pycisTopic.clust_vis import *
from pycisTopic.lda_models import *
import pyranges
import gzip
print('finished importing modules')

print()

print('functions...')
def make_argument_parser():
    """
    Creates an ArgumentParser to read the options for this script from
    sys.argv
    """
    parser = argparse.ArgumentParser(
        description=" Evaluate topic models.",)
    parser.add_argument('--input_cto', '-i', type=str, required=True,
                        help='Path to cisTopic object pickle file.')
    
    parser.add_argument('--topic_model_dir', '-tmd', type=str, required=True,
                        help='Path topic model dir.')

    parser.add_argument('--selected_model', '-m', type=int, required=True,
                        help='Path topic model dir.')

    parser.add_argument('--output_dir', '-o', type=str, required=True,
                        help='Path to output dir.')
    
    return parser


def load_and_add_topic_model(model_dir, model, cistopic_obj):
    """
    Load and add topic model
    """

    model_name = f'Topic{model}.pkl.gz'
    print(model_name)
    selected_model = pickle.load(gzip.open(os.path.join(model_dir, model_name),"rb"))
    cistopic_obj.add_LDA_model(selected_model)

def clustering_and_dim_red(cistopic_obj, model):
    """
    Finding clusters and running UMAP
    """
    ntopics = str(model)
    find_clusters(cistopic_obj,
              target  = 'cell',
              k = 10,
              res = [0.6, 1.2, 2.0, 2.5, 3.0, 3.5],
              prefix = f'{ntopics}_pycisTopic_',
              scale = True,
              split_pattern = '-')
    print('Finding clusters finished')
    print()
    
    print('Run UMAP.....')
    run_umap(cistopic_obj, target  = 'cell', scale=True, reduction_name = f'UMAP_{ntopics}')
    print('UMAP done!')
    
    print('Run tSNE.....')
    run_tsne(cistopic_obj, target  = 'cell', scale=True, reduction_name = f'tSNE_{ntopics}')
    print('tSNE done!')
    print()
    
    print('run harmony')
    harmony(cistopic_obj, vars_use = ['method'])
    print('harmony done')
    
    print('Run UMAP with harmony...')
    ## Run umap and tsne with harmony
    run_umap(cistopic_obj, target  = 'cell', scale=True, reduction_name = f'UMAP_{ntopics}_harmony', harmony = True)
    print('UMAP done!')
    
    print('run tSNE with harmony...')
    run_tsne(cistopic_obj, target  = 'cell', scale=True, reduction_name = f'tSNE_{ntopics}_harmony', harmony = True)
    print('tSNE done! ')

print()

print('starting main')
print()
def main():
    """
    The main executable function
    """
    print('make argument parser')
    parser = make_argument_parser()
    args = parser.parse_args()

    print('Opening CTO...')
    filename = args.input_cto
    infile = open(filename, 'rb')
    cistopic_obj = pickle.load(infile)
    infile.close()
    print('Input cisTopic_object: ', filename)
    print('Opening done')
    print()
    
    model_dir = args.topic_model_dir
    print('Model dir: ', model_dir)
    print()
    
    model = args.selected_model
    print('Selected topic model: ', str(model))
    print()
    
    output = args.output_dir
    print('Output dir:', output)
    if not os.path.exists(output):
        print('Output dir does not exist, creating...')
        os.makedirs(output)
        print(f"Directory '{output}' created.")
    else:
        print(f"Directory '{output}' already exists.")
    print()
   
    print('Load and select model:')
    load_and_add_topic_model(model_dir, model, cistopic_obj)
    print('model added')
    print()

    print('Clustering and UMAP....')
    clustering_and_dim_red(cistopic_obj, model)
    print()

    print('saving CTO...')
    ntopics = str(model)
    with open(output + f'/{ntopics}_cto.pkl', 'wb') as f:  
        pickle.dump(cistopic_obj, f)
    print('saved')
    print('___ FINISHED ___')


if __name__ == "__main__":
    main()