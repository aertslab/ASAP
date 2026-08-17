import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import polars as pl
import pickle as pk
from utils import *
import numpy as np
import argparse


def main():

    print()
    print(" NB dispersion parameter test for peak count in loci VS null-matched loci")
    print()
    
    varlink = "/lustre1/project/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/2026_full_dataset/20260703_caQTL_full/full_caqtl1kb_finemap200kb_30PCs/susie/credible_set_summary/credible_set_variants_combined.tsv.gz"
    variants = pl.read_csv(varlink, separator='\t')
    print(f" [peaks] Loaded {variants.shape[0]} variant scores from {varlink}")

    pip5 = variants.filter(pl.col('pip') > 0.5)
    print(" [peaks] Subsetted peaks on pip > 0.5")

    null_regions = pl.read_csv('6_enrichment/regions.tsv', separator='\t')
    lead_regions = pl.read_csv('4_regions/regions.tsv', separator='\t')

    with open ('6_enrichment/filtered_leads.pk', 'rb') as fp:
        filtered_leads = pk.load(fp)

    lead_regions_filtered = lead_regions.filter(pl.col("lead_id").is_in(filtered_leads))

    draw_dict = pk.load(open("6_enrichment/draw_dict.pkl","rb"))

    # Stratify peak count distributions per tissue & cell type
    unique_tissues = pip5[['region','cell_type']].unique()
    for tissue, cell_type in unique_tissues.filter(pl.col("cell_type").is_in(["Oligo","Micro-PVM"])).iter_rows():
        section = pip5.filter( (pl.col('cell_type') == cell_type) & (pl.col('region') == tissue) )
        peaks = parse_peaks(set(section['phenotype_id']))
        result = nb_dispersion_enrichment(peaks, lead_regions_filtered, null_regions, draw_dict)
        print(f" [peaks]  - {tissue}\t {cell_type}: identical NB distribution of peak counts in ref loci vs null loci 0-hypothesis p-value: {result['p_value']}")


if __name__ == "__main__":
    main()