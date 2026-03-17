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


# Ensure output directories exist
output_dirs = [
    'cc/ct/parquets/snp', 'cc/ct/results/snp',
    'cc/ct/parquets/indel', 'cc/ct/results/indel',
]
for d in output_dirs:
    os.makedirs(d, exist_ok=True)

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
def run_eqtl_analysis(plink_prefix_path, expression_bed, covariates_file, prefix, output_file):
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
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, prefix, covariates_df=covariates_df)

    # exmpirical p-values all genes
    cis_df = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df=covariates_df)


    cis_df.to_csv(output_file, index = True, compression='gzip')


# Specify the directory path
bed_path = '/cc/ct/phenotypes/'
suffix = '_log1p_norm.bed.gz' # suffix of input files


# Get a list of all files in the directory
all_files = [f for f in os.listdir(bed_path) if os.path.isfile(os.path.join(bed_path, f))]

print(all_files)

# define paths to data
# SNPs
plink_prefix_path = '/plink/WGS_chm13.snps_only'

# SNP eQTL analysis no PCA
for file in all_files:
    print(f'Now running for {file}')
    name = extract_name(file, suffix)
    print(name)
    expression_bed = os.path.join(bed_path, file)
    covariates_file = f'/cc/ct/covariates/snp/{name}_cov.txt'
    prefix = f'cc/ct/parquets/snp/{name}'
    output_file = f'cc/ct/results/snp/eQTL_{name}.csv.gz'
    print(f'Running for {name}')
    try:
        run_eqtl_analysis(plink_prefix_path, expression_bed, covariates_file, prefix, output_file)
    except Exception as e:
        print(f"An error occurred while running eQTL analysis for {name}: {e}")
    print("Done")




# time the script
end = datetime.datetime.now()
print("Time taken:", end-start)
print("All done!")
