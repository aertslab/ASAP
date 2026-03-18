## tensorQTL parquet files
print('importing...')
import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis, trans, post
import datetime
import argparse

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas: {pd.__version__}")
print(f"tensorqtl: {tensorqtl.__version__}")

print('importing done')
print()

## Arguments
def make_argument_parser():
    """
    Creates an ArgumentParser to read the options for this script from
    sys.argv
    """
    parser = argparse.ArgumentParser(
        description="Load assigned topics and save topic model",)
    parser.add_argument('--br', '-br', type=str, required=True,
                        help='brain region')
    parser.add_argument('--ct', '-ct', type=str, required=True,
                        help='cell type')
    parser.add_argument('--vt','-vt', type=str, required=True, 
                        help='variant type (SNP/Indel/SV)')              
    parser.add_argument('--cov','-cov', type=str, required=True, 
                        help='covariates to use (no_pca/genotype_pca/genotype_data_pca)')
    parser.add_argument('--peaks','-peaks', type=str, required=True, 
                        help='Peak counts to use (CT1/CT2)')
    return parser


def get_file_paths(ct, vt, cov, peaks):

    if vt == 'SNP':
        plink_prefix_path = '/data/plink/WGS_chm13.snps_only'
    elif vt == 'Indel':
        plink_prefix_path = '/data/plink/WGS_chm13.indels_only'
    elif vt == 'SV':
        plink_prefix_path = '/data/plink/WGS_chm13_ASAP_SV'
    else:
        print('wrong variant type argument')

    
    if peaks == 'CT1':
        break
    elif peaks == 'CT2':
        expression_bed = f'/snATAC_analysis/{br}/caQTL/data/{ct}_counts_summed_per_donor_sorted_cpm_log1p.bed'
        if cov == 'no_pca':
            covariates_file = f'/snATAC_analysis/{br}/caQTL/data/covariates/{ct}_cov.txt'
        elif cov == 'genotype_pca':
            covariates_file = f'/snATAC_analysis/{br}/caQTL/data/{ct}_cov_{vt}.txt'
        elif cov == 'genotype_data_pca_30':
            covariates_file = f'/snATAC_analysis/{br}/caQTL/data/{ct}_cov_{vt}_data_pca_30PCs.txt'
        else:
            print('wrong cov argument')
    else:
        print('wrong peaks argument')

    prefix = f'/snATAC_analysis/{br}/caQTL/out/parquets/{peaks}_peak_regions/{cov}/{ct}_{vt}' 

    print(f'plink: {plink_prefix_path}')
    print(f'expression_bed: {expression_bed}')
    print(f'covariates_file: {covariates_file}')
    print(f'prefix: {prefix}')
    print()
    
    return plink_prefix_path, expression_bed, covariates_file, prefix

def make_parquet(ct, vt, cov, peaks):

    if vt == 'SNP':
        w = 0
    elif vt == 'Indel':
        w =  200   
    elif vt == 'SV':
        w = 3000
    else:
        print('wrong variant type argument') 

    plink_prefix_path, expression_bed, covariates_file, prefix = get_file_paths(ct, vt, cov, peaks)

    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
    
    # PLINK reader for genotypes
    pgr = pgen.PgenReader(plink_prefix_path)
    genotype_df = pgr.load_genotypes()
    variant_df = pgr.variant_df
    
    columns_to_keep = [col for col in phenotype_df.columns if col in genotype_df.columns]
    phenotype_df = phenotype_df[columns_to_keep]
    covariates_df = covariates_df.loc[columns_to_keep]
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, covariates_df=covariates_df, window=w)

def main():
    """
    The main executable function
    """
    print('make argument parser')
    parser = make_argument_parser()
    args = parser.parse_args()
    
    ## Print the arguments
    ct = args.ct
    print('Cell type: ', ct)

    vt = args.vt
    print('variant type: ', vt)

    cov = args.cov
    print('covariates: ', cov)

    peaks = args.peaks
    print('peaks: ', peaks)

    print('get parquet files...')
    make_parquet(ct, vt, cov, peaks)
    print('done')
    print('-------------------------')
  

if __name__ == "__main__":
    main()    