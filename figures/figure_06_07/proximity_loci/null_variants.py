from utils import match_distribution
import pandas as pd
import argparse

def main():

    args = parse_args()

    leads = pd.read_csv('2_source/GP2_leads.tsv', sep='\t')
    leads['CHR'] = leads['CHR'].astype(int)
    leads['AF'] = leads['AF'].astype(float)
    leads['MAF'] = leads.apply(lambda r: r['AF'] if r['AF'] < 0.5 else 1 - r['AF'], axis=1)

    variants = pd.read_csv('1_plink/1KG.T2T.tsv', sep='\t')
    variants = variants[variants['MAF'] > 0.01]
    variants['CHR'] = [int(chr.replace('chr','')) for chr in variants['CHR']]

    null_variants = match_distribution(leads, variants, chrom_col="CHR", maf_col="MAF", n_bins=10, oversample_factor=args.sample)
    null_variants['BP'] = null_variants.apply(lambda r: r['ID'].split('_')[1], axis=1)
    null_variants['Nearest Gene'] = 'Irrelevant'
    null_variants[['CHR','BP','Nearest Gene','ID','MAF']].to_csv('6_enrichment/null_variants.tsv', sep='\t', index=None)


def parse_args():
    parser = argparse.ArgumentParser(
        description="",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--sample", type=int, default=100,
                         help="Sample 'sample' x len(lead-variants) variants")
    return parser.parse_args()

if __name__ == "__main__":
    main()