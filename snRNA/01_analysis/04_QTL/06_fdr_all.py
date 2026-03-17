import pandas as pd
import polars as pl

import io
import sys
import glob
import re

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import os

from statsmodels.stats.multitest import multipletests
import statsmodels

# cell type name correction
cell_type_labels_corr={"L2_3_IT": "L2/3 IT","L5_IT": "L5 IT", "L6_CT": "L6 CT", "L4_IT": "L4 IT", "L5_IT_A": "L5 IT A", "L5_IT_B": "L5 IT B", "L6_IT": "L6 IT", "L6_IT_Car3": "L6 IT Car3", "L5_6_NP": "L5/6 NP", "neurons" : "Neurons"}



def process_directory(path, directory, prefix_list):
    cis_dfs = []

    for prefix in prefix_list:
        prefix_res = f'eQTL_{prefix}'
        
        chr_files = {
        os.path.basename(f).split('.')[2]: (
            pl.read_parquet(f)
        )
            for f in glob.glob(f'{path}/parquets/{directory}/{prefix}.cis_qtl_pairs.*.parquet')
        }
        
        cell_type = prefix.split('_')[2:]
        cell_type = "_".join(cell_type)
        # Concatenate Polars DataFrames
        cis_df = pl.concat(list(chr_files.values()))
        cis_df = cis_df.with_columns([
            pl.lit(cell_type).alias("cell_type"),
            pl.lit(directory).alias("variant_type")
        ])
        cis_dfs.append(cis_df)

    return cis_dfs

def process_all_directories(path, directories, prefix_list):
    all_cis_dfs = []

    for directory in directories:
        cis_dfs = process_directory(path, directory, prefix_list)
        all_cis_dfs.extend(cis_dfs)

    # Concatenate all Polars DataFrames
    combined_cis_df = pl.concat(all_cis_dfs)
    return combined_cis_df
    

directories = ['snp', 'indel']

prefix_list = prefix_list = ['20260219_sn_Astro',
    '20260219_sn_Micro-PVM',
    '20260219_sn_OPC',
    '20260219_sn_Oligo',
    '20260219_sn_neurons',]


path = '/sn/'
# Process all directories
combined_cis_df_sn = process_all_directories(path, directories, prefix_list)
combined_cis_df_sn = combined_cis_df_sn.filter(~pl.col("pval_nominal").is_null())

combined_cis_df_sn = combined_cis_df_sn.with_columns(
    pl.col("cell_type").replace(cell_type_labels_corr).alias("cell_type")
)

# Group by cell_type and variant_type, count unique phenotype_id
result = (
    combined_cis_df_sn
    .group_by(["cell_type", "variant_type"])
    .agg(pl.col("phenotype_id").n_unique().alias("phenotype_id_nunique"))
)

print(result)


## CC
    


directories = ['snp', 'indel']

prefix_list = prefix_list = [
    '20260219_cc_Astro',
'20260219_cc_L2_3_IT',
'20260219_cc_L4_IT',
'20260219_cc_L5_6_NP',
'20260219_cc_L5_IT_A',
'20260219_cc_L5_IT_B',
'20260219_cc_L6_CT',
'20260219_cc_L6_IT',
'20260219_cc_L6b',
'20260219_cc_Lamp5',
'20260219_cc_Micro-PVM',
'20260219_cc_OPC',
'20260219_cc_Oligo',
'20260219_cc_Pvalb',
'20260219_cc_Sncg',
'20260219_cc_Sst',
'20260219_cc_Vip',
]



path = '/cc/ct/'
# Process all directories
combined_cis_df_cc = process_all_directories(path, directories, prefix_list)
combined_cis_df_cc = combined_cis_df_cc.filter(~pl.col("pval_nominal").is_null())

combined_cis_df_cc = combined_cis_df_cc.with_columns(
    pl.col("cell_type").replace(cell_type_labels_corr).alias("cell_type")
)

# Group by cell_type and variant_type, count unique phenotype_id
result = (
    combined_cis_df_cc
    .group_by(["cell_type", "variant_type"])
    .agg(pl.col("phenotype_id").n_unique().alias("phenotype_id_nunique"))
)

print(result)


## All variants together
combined_cis_df_cc = combined_cis_df_cc.with_columns(
    pl.lit("CC").alias("region")
)

combined_cis_df_sn = combined_cis_df_sn.with_columns(
    pl.lit("SN").alias("region")
)
combined_cis_df = pl.concat([combined_cis_df_cc, combined_cis_df_sn])
combined_cis_df = combined_cis_df.with_columns(
    (pl.col("cell_type") + ", " + pl.col("region")).alias("ct_reg")
)



combined_cis_df = combined_cis_df.filter(combined_cis_df['pval_nominal'].is_not_null())


p_values = combined_cis_df["pval_nominal"].to_list()

# Perform FDR correction
signif, p_adj = statsmodels.stats.multitest.fdrcorrection(p_values, alpha=0.05, is_sorted=False)


combined_cis_df = combined_cis_df.with_columns([
    pl.Series("signif", signif),  
    pl.Series("p_adj", p_adj)    
])


combined_cis_df.write_parquet("out/eQTL_2026_protein_coding_all.parquet")