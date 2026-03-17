import decoupler as dc
import pertpy as pt
import scanpy as sc
from natsort import natsorted
import pandas as pd
import numpy as np
import gzip
import datetime
from sklearn.decomposition import PCA
from sklearn import preprocessing
import os


import warnings
warnings.filterwarnings("ignore")


date = datetime.datetime.now().strftime('%Y%m%d')

start = datetime.datetime.now()

# make sure output dirs exist
output_dirs = [
    'cc/ct/phenotypes',
    'cc/ct/covariates/snp',
    'cc/ct/covariates/indel',
    'cc/ct/covariates'
]

for d in output_dirs:
    os.makedirs(d, exist_ok=True)

print("Today is",date)
print("Starting the script")



#####################################
# define function to load gtf files
#####################################
def process_gtf(file_path):
    # Initialize lists to store the data
    records = []
    
    # Open and read the gzipped file
    with gzip.open(file_path, 'rt') as f:
        for line in f:
            if line.startswith('#'):  # Skip header lines
                continue
            
            fields = line.strip().split('\t')
            
            # Check if third column is 'gene'
            if fields[2] == 'gene':
                chr_num = fields[0]
                start = fields[3]
                end = fields[4]
                
                # Process column 9 to extract gene ID
                attributes = fields[8]
                gene_id = attributes.split(';')[0]  # Get first attribute
                gene_id = gene_id.split('"')[1]  # Extract value between quotes
                gene_id = gene_id.replace("_", "-")  # IMPORTANT since cellranger seem to do this internally which led to some gene drop out

                # process variable field in columns 9 to get gene biotype and keep only protein coding
                # of note, MT genes are also excluded because they have gene_type instead of gene_biotype
                cells = attributes.split(';')
                # Initialize biotype
                biotype = None

                for cell in cells:
                    if 'gene_biotype' in cell:
                        biotype = cell.split('"')[1]
                        break

                # Check if biotype is found and print it
                if biotype:
                    #print(biotype)
                    if biotype == 'protein_coding':
                        records.append([chr_num, start, end, gene_id])
    
    # Create DataFrame
    df = pd.DataFrame(records, columns=['#chr', 'start', 'end', 'gene_id'])
    return df



################################
# Function for PCA calculation from expression data
# (for covariates)
################################

def perform_pca(pdata, n_components=30, layer='log1p', cell_type='Micro-PVM'):
    """
    Perform PCA on the given pdata object and return the resulting PCA DataFrame.

    Parameters:
    pdata (AnnData): The AnnData object containing the data.
    n_components (int): The number of PCA components to compute.
    layer (str): The layer of the AnnData object to use for PCA.
    cell_type (str): The cell type to be replaced in the column names.

    Returns:
    pd.DataFrame: The resulting PCA DataFrame.
    """
    # get data from wanted layer
    counts_df = pdata.to_df(layer=layer)
    
    # PCA
    pca = PCA(n_components=n_components)
    pca.fit(counts_df)
    pca_data = pca.transform(counts_df)
    
    # results --> add expr to name
    # results --> transpose to fit in covariate space
    pca_df = pd.DataFrame(pca_data, index=counts_df.index, columns=[f'PC{i+1}_expr' for i in range(n_components)])
    pca_df.index = pca_df.index.str.replace(f"_{cell_type}", "", regex=False)
    pca_df = pca_df.T
    
    return pca_df


################

# load gtf file
genes = process_gtf('/staging/leuven/res_00001/genomes/homo_sapiens/CHM13v2_maskedY_rCRS/indexes/cellranger/8.0.1/CHM13v2_maskedY_rCRS/genes/genes.gtf.gz')

# load PCA file for SNPs
snps_donor_pca = '../plink/WGS_chm13.snps_only_genotype_pca.eigenvec'
snps_donor_pca_df = pd.read_csv(snps_donor_pca, delim_whitespace=True)
snps_donor_pca_df.index = snps_donor_pca_df['#IID']
snps_donor_pca_df = snps_donor_pca_df.drop(columns=['#IID'])
snps_donor_pca_df = snps_donor_pca_df.T
snps_donor_pca_df.index = [f'SNP_{x}' for x in snps_donor_pca_df.index]

# load PCA file for INDELS
indel_donor_pca = '../plink/WGS_chm13.indels_only_genotype_pca.eigenvec'
indel_donor_pca_df = pd.read_csv(indel_donor_pca, delim_whitespace=True)
indel_donor_pca_df.index =indel_donor_pca_df['#IID']
indel_donor_pca_df = indel_donor_pca_df.drop(columns=['#IID'])
indel_donor_pca_df = indel_donor_pca_df.T
indel_donor_pca_df.index = [f'INDEL_{x}' for x in indel_donor_pca_df.index]

# load adata
path_to_adata = '../02_combination_10x_pb/hdf5/20250811_all_data_integrated_filtered_umap_scANVI_metadata_ct_corrected.annot_corrected.index_dedup.cg.h5ad'
adata = sc.read_h5ad(path_to_adata)
adata.obs.set_index(adata.obs.new_index.astype(str),  inplace=True)

# prepare data for pseudobulk
adata.X = adata.layers["counts"]
adata = adata[~adata.obs['experiment'].isin(['mo005', 'mo009'])]
adata = adata[~adata.obs['region'].isin(['Motor Cortex'])]
adata.obs.loc[adata.obs['primary_diagnosis'] == 'Other neurological disorder', 'primary_diagnosis'] = 'Healthy Control'
adata.obs['primary_diagnosis'] = adata.obs['primary_diagnosis'].cat.rename_categories({'Healthy Control': 'Control', 'Idiopathic PD': 'PD'})

# add info on different batches
# some donors have data from both PB and 10X and we have to take that into acccount for modelling
grouped = adata.obs.groupby(['vireo_assignment', 'batch']).size().reset_index(name='counts')
combined_dict = {}
for donor, group in grouped.groupby('vireo_assignment'):
    if group['counts'].gt(0).sum() > 1:
        combined_dict[donor] = 'mix'
    else:
        combined_dict[donor] = group.loc[group['counts'].gt(0), 'batch'].values[0]
adata.obs['batch_combined'] = adata.obs['vireo_assignment'].map(combined_dict)

# create pseudobulk
ps = pt.tl.PseudobulkSpace()
pdata = ps.compute(adata, target_col="vireo_assignment", groups_col="consensus_cell_type_corrected", layer_key="counts", mode="sum", min_cells=50, min_counts=1000)

# get cell type specific pseudobulk
cell_types = adata.obs['consensus_cell_type_corrected'].cat.categories

# print number of unique cell types
print(f"Number of unique cell types: {len(cell_types)}")

for cell_type in cell_types:
    print(f"Processing {cell_type} pseudobulk expression data") 
    print(f"Number of unique donors: {pdata[pdata.obs['consensus_cell_type_corrected'] == cell_type].obs['vireo_assignment'].nunique()}")

    if pdata[pdata.obs['consensus_cell_type_corrected'] == cell_type].obs['vireo_assignment'].nunique() < 30:
        print(f"Skipping {cell_type} because of low number of donors")
        continue
   
    cell_type_data = pdata[pdata.obs['consensus_cell_type_corrected'] == cell_type]
    cell_type_data.layers["counts"] = cell_type_data.X.copy()
    sc.pp.normalize_total(cell_type_data, target_sum=1e6) # target increaserc so it is CPM norm
    sc.pp.log1p(cell_type_data)
    cell_type_data.layers["log1p"] = cell_type_data.X.copy()

    # save the log normalised expression data
    log1p_df = cell_type_data.to_df(layer='log1p')
    log1p_df = log1p_df.T
    log1p_df.columns = log1p_df.columns.str.replace(f'_{cell_type}', "", regex=False)
    log1p_df.reset_index(inplace=True)
    log1p_df.rename(columns={'index': 'gene_id'}, inplace=True)

    # merge expression with gene information
    final_df = pd.merge(genes, 
         log1p_df, 
         on='gene_id',
         how='inner')
    
    print(f"Number of genes: {final_df.shape[0]}")

    # prepare output
    final_df['start'] = pd.to_numeric(final_df['start'], errors='coerce')
    final_df['end'] = pd.to_numeric(final_df['end'], errors='coerce')
    final_df['#chr'] = pd.Categorical(final_df['#chr'], categories=natsorted(final_df['#chr'].unique()), ordered=True)

    final_df = final_df.sort_values(by=['#chr', 'start','end'])

    # if cell type name has / replace with _
    cell_type_out = cell_type.replace("/", "_")
    cell_type_out = cell_type_out.replace(" ", "_")
    output_path = f'cc/ct/phenotypes/{date}_cc_{cell_type_out}_log1p_norm.bed.gz'
    final_df.to_csv(output_path,index = False, sep = "\t", compression='gzip')

    print(f"Saved {cell_type} pseudobulk expression data to {output_path}")

    # get covariates
    cov  = cell_type_data.obs
    columns_to_keep = ['biobank_name', 'sex','age_at_collection','batch_combined']
    # Select columns to keep using loc
    cov = cov.loc[:, columns_to_keep]
    cov = pd.get_dummies(cov, columns= ['biobank_name','sex','batch_combined'], drop_first=True)
    cov = cov.astype(int)
    cov = cov.T
    cov.columns = cov.columns.str.replace(f'_{cell_type}', "", regex=False)
    


    cov_snp = pd.concat([cov, snps_donor_pca_df], axis=0)
    cov_snp.to_csv(f'cc/ct/covariates/snp/{date}_cc_{cell_type_out}_cov.txt',sep='\t')

    cov_indel = pd.concat([cov, indel_donor_pca_df], axis=0)
    cov_indel.to_csv(f'cc/ct/covariates/indel/{date}_cc_{cell_type_out}_cov.txt',sep='\t')


    # add expression pca
    # get PCA covariates
    pca_df = perform_pca(cell_type_data, n_components=30, layer='log1p', cell_type=cell_type)
    # combine covariates
    cov = pd.concat([cov, pca_df], axis=0)

    cov_snp = pd.concat([cov, snps_donor_pca_df], axis=0)
    cov_snp.to_csv(f'cc/ct/covariates/snp/{date}_cc_{cell_type_out}_pca_cov.txt',sep='\t')

    cov_indel = pd.concat([cov, indel_donor_pca_df], axis=0)
    cov_indel.to_csv(f'cc/ct/covariates/indel/{date}_cc_{cell_type_out}_pca_cov.txt',sep='\t')

    # for plotting - add PD as well 
    cov  = cell_type_data.obs
    columns_to_keep = ['biobank_name', 'sex','age_at_collection','primary_diagnosis','batch_combined']
    # Select columns to keep using loc
    cov = cov.loc[:, columns_to_keep]
    cov = pd.get_dummies(cov, columns= ['biobank_name','sex','primary_diagnosis','batch_combined'], drop_first=True)
    cov = cov.astype(int)
    cov = cov.T
    cov.columns = cov.columns.str.replace(f'_{cell_type}', "", regex=False)
    cov.to_csv(f'cc/ct/covariates/{date}_cc_{cell_type_out}_noexpr_cov_plotting.txt',sep='\t')
    print(f"Saved {cell_type} covariates to covariates/{date}_{cell_type_out}_cc_cov.txt")

    # cell type done, move to next
    print(f"Done with {cell_type} pseudobulk expression data")
    print("")
    print("")



# time the script
end = datetime.datetime.now()
print("Time taken:", end-start)
print("All done!")



    
