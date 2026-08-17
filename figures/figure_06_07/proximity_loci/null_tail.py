from scipy.stats import mannwhitneyu
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
    print(" Tail-fatness test on CREsted scores in loci vs null-matched loci")
    print()

    args = parse_args()

    if '.tsv' in args.score_dataframe:
        variants = pl.read_csv(args.score_dataframe, separator='\t')
    elif '.csv' in args.score_dataframe:
        variants = pl.read_csv(args.score_dataframe)
    elif '.parquet' in args.score_dataframe:
        variants = pl.read_parquet(args.score_dataframe)
    else:
        print(f" [tail] Error: Implement reading this filetype: {args.score_dataframe}")
        return
    print(f" [tail] Loaded {variants.shape[0]} variant scores from {args.score_dataframe}")

    null_regions = pl.read_csv('6_enrichment/regions.tsv', separator='\t')
    lead_regions = pl.read_csv('4_regions/regions.tsv', separator='\t')

    with open ('6_enrichment/filtered_leads.pk', 'rb') as fp:
        filtered_leads = pk.load(fp)

    lead_regions_filtered = lead_regions.filter(pl.col("lead_id").is_in(filtered_leads))

    draw_dict = pk.load(open("6_enrichment/draw_dict.pkl","rb"))

    # score distribution of variants in lead blocks pooled VS score distribution of variants in null loci pooled
    stat, p = mwu_lead_vs_pooled_null_variants(variants.sample(min([variants.shape[0], 1_000_000])), lead_regions_filtered, null_regions, draw_dict, args.metric_column)
    print(f" [tail] Variants in lead loci distribution identical to pooled null loci leads 0-hypothesis: p: {p}")

    # Take abs value and cutof variant scores
    variants_filtered = variants.with_columns(pl.col([args.metric_column]).abs().alias('metric'))
    variants_filtered = variants_filtered.filter(pl.col('metric') > args.threshold)
    print(f" [tail] {variants_filtered.shape[0]} variants kept after filtering {args.metric_column} > {args.threshold}")

    result = gpd_tail_enrichment(variants_filtered, lead_regions_filtered, null_regions, draw_dict,
                                   value_col=args.metric_column, threshold=args.threshold)
    print(f" [tail] Tail fatness of lead loci variant {args.metric_column} is fatter than those of null matched loci sets 0-hypothesis: p: {result['p_value']}")



def parse_args():
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--score-dataframe", type=str, required=True,
                         help="tsv of variants score according to a certain metric")
    parser.add_argument("--metric-column", type=str, required=True,
                         help="column name of metric of score_dataframe")
    parser.add_argument("--threshold", type=float, required=True,
                         help="Threshold of the value above which the metric shoulkd be kept")
    return parser.parse_args()

if __name__ == "__main__":
    main()