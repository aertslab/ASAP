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
    print(" Calculating loci parameters and bootstrapping loci sets to match lead-loci")
    print()

    args = parse_args()

    # Load lead and null loci 1KG SNP overlap
    lead_snps = pd.read_csv('6_enrichment/lead_loci_SNP.tsv', sep='\t')
    null_snps = pd.read_csv('6_enrichment/null_loci_SNP.tsv', sep='\t')

    # calculate metrics (SNP density, mean MAF, size...) per locus
    lead_metrics_raw = build_locus_metrics(lead_snps)
    lead_metrics_raw['SNP_density'] = lead_metrics_raw['snp_count'] / lead_metrics_raw['size']

    # Filter lead loci based on (hardcoded) outlier thesholds
    filtered_leads = set(lead_metrics_raw[(lead_metrics_raw['mean_maf'] < 0.1) & (lead_metrics_raw['snp_count'] < 80_000) & (lead_metrics_raw['size'] < 2_500_000)]['locus_id'])
    lead_snps_filtered = lead_snps[lead_snps['locus_id'].isin(filtered_leads)]
    lead_snps_filtered = lead_snps_filtered[~lead_snps_filtered['locus_id'].isin(['chr6_70943119','chr15_86728586','chr21_38965849','chr10_120548988'])]

    # Save list of non filtered leads
    with open('6_enrichment/filtered_leads.pk', 'wb') as fp:
        pk.dump(set(lead_snps_filtered['locus_id']), fp)
    print(f" [loci] Saved lead loci sets to 6_enrichment/filtered_leads.pkl after filtering ({len(set(lead_snps_filtered['locus_id']))} remainig)")

    # Load peak df
    peaks = pd.read_csv(args.peaks_df, sep='\t', header=None, names=['chr','start','end','name','height'])

    # Calculate metrics again after filtering, also for null loci
    lead_metrics = build_locus_metrics(lead_snps_filtered, peaks_df=peaks)
    null_candidate_metrics = build_locus_metrics(null_snps, peaks_df=peaks)

    # Sample null loci to match stratified lead loci distribution
    null_draws, shortfalls, bin_edges = sample_null_regions(
            lead_metrics, null_candidate_metrics, 
            linear_cols=("mean_maf","SNP_density","peak_density"), n_bins=3, M=100, random_state=0)

    # Plot metric densities for visual check, calculate distributions difference p_value
    for metric in ['mean_maf','SNP_density','size','peak_density']:
        stat, p = mwu_lead_vs_pooled_null(lead_metrics, null_candidate_metrics, null_draws, col=metric)
        print(f" [loci] Mann-Withney-U test for lead loci metric distribution being equal to pooled null loci metric distribution of {metric} 0-hypothesis: p={p}")
        plot_strat_density(lead_metrics, null_candidate_metrics, null_draws, metric, p_val=p, bin_edges=bin_edges)

    # Save bootstraps of null loci
    draw_dict = dict(zip(range(len(null_draws)),null_draws))
    pk.dump(draw_dict, open("6_enrichment/draw_dict.pkl","wb"))
    print(f" [loci] Saved bootstrap sets to 6_enrichment/draw_dict.pkl")


def parse_args():
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--peaks-df", type=str, required=True,
                         help="bed file of consensus peaks to match loci (chr, start, end, name, height)")
    return parser.parse_args()


if __name__ == "__main__":
    main()