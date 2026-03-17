import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pyranges as pr

import sys
import argparse


# script to calculate average TSS methylation for top and bottom 10% expressed genes in a given cell type, based on a bedmethyl file and a TSS bed file.

# Check if the correct number of arguments is provided
parser = argparse.ArgumentParser(description="Calculate TE methylation from bedmethyl file.")
parser.add_argument("--sample", required=True, help="Path to the bedmethyl file")
parser.add_argument("--cell_type", default="oligo", help="Cell type to filter the bedmethyl file (default: oligo)")
args = parser.parse_args()

bedmethyl_path = args.sample
cell_type = args.cell_type


# check if file exists
if not os.path.isfile(bedmethyl_path):
    print(f"File {bedmethyl_path} does not exist.")
    sys.exit(1)


# define the function to calculate TE methylation
def tss_methylation(bedmethyl_path, tss, cell_type = "oligo"):

    # load in methylation data
    methdata = pd.read_csv(
        bedmethyl_path, sep='\t',
        header=None,
        index_col=False,
        names=["Chromosome", "Start", "End", "name", "score", "strand", "tstart", "tend", "color", "coverage", "freq", "mod", "canon","other","delete","fail","diff","nocall"])

    # keep only m
    methdata = methdata[methdata.name == 'm']
    # coverage over 5 
    methdata = methdata[methdata.coverage > 5]

    # add info on expression - so we only check TSS at Top and bottom genes
    expression_low = pd.read_csv(f'expression_info/cc_{cell_type}_bottom10_expressed_genes.csv')  # oligo...
    expression_low['expression'] = 'low'
    expression_high = pd.read_csv(f'expression_info/cc_{cell_type}_top10_expressed_genes.csv')
    expression_high['expression'] = 'high'
    expression = pd.concat([expression_high, expression_low], ignore_index=True)
    
    tss = pd.merge(tss, expression, how = 'left', right_on = 'gene', left_on = 'gene_name')
    
    
    # convert to pyranges and overlap
    tss_r = pr.PyRanges(tss)
    tss_r = tss_r.slack(2000)
    methdata_r = pr.PyRanges(methdata)
    # overlap the two
    overlap = methdata_r.join(tss_r, how = None)
    overlap_df = overlap.df # and back to pandas

    overlap_df_filtered = overlap_df.dropna(subset=['expression']) # get rid of TSS not in top and bottom 10%
    result = overlap_df_filtered.groupby('expression')['freq'].mean().reset_index()

    result.columns = ['expression', 'mean_methylation']
    
    return result





tss_path = "/staging/leuven/res_00001/genomes/homo_sapiens/CHM13v2_maskedY_rCRS/annotations/CHM13v2_maskedY_rCRS.tss.bed"  # change to the TSS path
# load in TE file 
tss = pd.read_csv(
        tss_path, sep='\t', header=None, index_col=False,
        names=["Chromosome", "Start", "End", "gene_name", "x", "strand"])


filename = os.path.basename(bedmethyl_path)
# remove everythiing after . 

filename = filename.split('.')[0]

print(f'Now processing {filename}')
result = tss_methylation(bedmethyl_path, tss, cell_type)
result['sample'] = filename

print(f'{filename} done')

result.to_csv(f'results/{filename}_{cell_type}_methylation_tss.csv', index=False)
