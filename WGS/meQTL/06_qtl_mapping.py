import pandas as pd
import torch
import tensorqtl
from tensorqtl import pgen, cis, trans, post
import datetime
import os
import gzip


start = datetime.datetime.now()
date = datetime.datetime.now().strftime('%Y%m%d')

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas: {pd.__version__}")
print(f"tensorqtl: {tensorqtl.__version__}")



def extract_name(file_path, suffix):
    """
    Extracts the desired part of the file name by removing the specified suffix.
    
    Parameters:
    file_path (str): The full path to the file.
    suffix (str): The suffix to remove from the file name.
    
    Returns:
    str: The extracted file name.
    """
    base_name = os.path.basename(file_path)
    if base_name.endswith(suffix):
        return base_name[:-len(suffix)]
    else:
        return base_name



# define function to run tensorQTL
def run_qtl_analysis(plink_prefix_path, expression_bed, covariates_file, prefix, output_file, window = 100000):
    # load phenotypes and covariates
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

    # PLINK reader for genotypes
    pgr = pgen.PgenReader(plink_prefix_path)
    genotype_df = pgr.load_genotypes()
    variant_df = pgr.variant_df

    # this is for SV - get rid of gene expression in donors for which there is no genotype available
    columns_to_keep = [col for col in phenotype_df.columns if col in genotype_df.columns]
    phenotype_df = phenotype_df[columns_to_keep]
    covariates_df = covariates_df.loc[columns_to_keep]

    # map all cis-associations (results for each chromosome are written to file)
    # nominal p-values all genes
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, covariates_df=covariates_df, window = window)

    # exmpirical p-values all genes
    cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df=covariates_df, window = window)


    cis_df.to_csv(output_file, index = True, compression='gzip')




# define paths to data
# SNPs
expression_bed = '/tensor_input_files/20251021_mQTL_atac_v1_phenotypes.noindex.bed.gz'

plink_prefix_path = '/plink/WGS_chm13.snps_only.sample_subset'
covariates_file = f'/tensor_input_files/mQTL_SNP_cov.txt'
# get date for output files


prefix = f'out/{date}_mqtl/{date}_mqtl_snp'
#check prefix exists
os.makedirs(os.path.dirname(f'out/{date}_mqtl/'), exist_ok=True)
output_file = f'out/{date}_mqtl/mQTL_results_snp.csv.gz'

run_qtl_analysis(plink_prefix_path, expression_bed, covariates_file, prefix, output_file, window = 10000)


#####
#
#
#indels

plink_prefix_path = '/plink/WGS_chm13.indels_only'
covariates_file = f'/tensor_input_files/mQTL_indel_cov.txt'

prefix = f'out/{date}_mqtl/{date}_mqtl_indels'
#check prefix exists
os.makedirs(os.path.dirname(f'out/{date}_mqtl/'), exist_ok=True)
output_file = f'out/{date}_mqtl/mQTL_results_indels.csv.gz'

run_qtl_analysis(plink_prefix_path, expression_bed, covariates_file, prefix, output_file, window = 10000)

    
# time the script
end = datetime.datetime.now()
print("Time taken:", end-start)
print("All done!")
