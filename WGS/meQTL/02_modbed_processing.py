#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyranges as pr
import argparse
import os


parser = argparse.ArgumentParser(description='Process some input arguments.')
parser.add_argument('--input_file', type=str, help='Path to the input .bed.gz file')

args = parser.parse_args()
input_file = args.input_file

# Extract sample name from path
basename = os.path.basename(input_file)
sample = basename.split('_cpg.bed.gz')[0]

# check if the file exists
if not os.path.isfile(input_file):
    print(f"File {input_file} does not exist.")
    exit(1)

# load genome binning regions
chrom_bins = pd.read_csv('CC_CT2_consensus_regions_genome_binning.tsv', sep="\t") # potentially update the path to your genome binning file
chrom_bins = pr.PyRanges(chrom_bins)

# load modbed file
modbed = pd.read_csv(
    input_file,
    sep='\t', header=None, compression='gzip'
)


modbed = modbed[[0, 1, 2, 3, 10]] # keep only relevant columns
modbed.columns = ['Chromosome', 'Start', 'End', 'Mod', 'Fraction_modified']
modbed = modbed[modbed['Mod'].str.contains('m')]  # drop hydroxymethylation rows
modbed = pr.PyRanges(modbed)

overlap = modbed.join(chrom_bins)
df = overlap.df.drop(columns=['Start', 'End','Mod'])

region_counts = df.groupby("region_id").size()
valid_regions = region_counts[region_counts >= 5].index

# Filter df to only valid region_ids - the ones that have at least 5 C with information
filtered_df = df[df["region_id"].isin(valid_regions)]

# combine + average the fraction modified within bin
fs = {"Chromosome": "first", "Start_b":
      "first", "End_b": "first", "Fraction_modified": "mean"}
result = filtered_df.groupby("region_id").agg(fs).reset_index()

result = result.rename(
    columns={"Chromosome": "Chromosome", "Start_b": "Start",
             "End_b": "End"})
# reorder 
column_order = ["Chromosome", "Start", "End", "Fraction_modified", "region_id"]
result = result[column_order]
result_r = pr.PyRanges(result)
result_r = result_r.sort()

output_dir = '/donor_binned_data' # potentially update the path to your output directory
os.makedirs(output_dir, exist_ok=True)
result_r.to_csv(f'{output_dir}/{sample}.bed.gz')
