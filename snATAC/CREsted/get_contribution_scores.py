import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import crested
#import tangermeme
import keras
import polars as pl
import re
import math
from crested.utils import one_hot_encode_sequence
import math
import numpy as np
from functools import lru_cache
analysis_dir = "/staging/leuven/stg_00090/ASA/analysis/"
model_dir = "/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/peak_regression/"
cell_types_corr={"L23_IT": "L2_3_IT", "L56_NP": "L5_6_NP"}


# run with:
# source /staging/leuven/stg_00002/lcb/osiga/software/anaconda3/etc/profile.d/conda.sh
# conda activate /lustre1/project/stg_00002/mambaforge/vsc35059/envs/crested_070525


### Functions ###


def extend_coord(start, end, extend_length=2114):
    """
    Extend the given sequence to a specific length.

    Args:
        start (int): Start position of the sequence.
        end (int): End position of the sequence.
        extend_length (int, optional): Length to extend the sequence.

    Returns:
        tuple: Updated start and end positions
    """
    diff = extend_length - (end - start)
    start = start - int(math.ceil(diff / 2))
    end = end + int(math.floor(diff / 2))

    return start, end


def split_coordinates_peak(region:str, extend_peak = False, extend_length=2114):
    
    match = re.search(r'(\w+):(\d+)-(\d+)', region)
    chrom, start, end = match.group(1), int(match.group(2)), int(match.group(3))
    
    if extend_peak:
        start, end = extend_coord(start, end, extend_length=extend_length)
    
    return (chrom, start, end)


def split_coordinates_variant(var:str):

    var_info = var_id.split("_")
    return (var_info[0], int(var_info[1]), var_info[2], var_info[3])


@lru_cache(maxsize=100000)
def cached_fetch(genome, chrom, start, end):
    return genome.fetch(chrom, start, end).upper()
    

def get_var_sequence(
    var_id,
    peak_id,
    genome,
    peak_width = 2114,
    check_ref_sequence=True
):

    chrom, var_pos, ref, alt = split_coordinates_variant(var_id)
    ref = ref.upper()
    alt = alt.upper()
    chrom2, peak_start, peak_end = split_coordinates_peak(peak_id, extend_peak = True, extend_length=peak_width)
    rel_pos = var_pos - peak_start - 1 
    assert chrom == chrom2

    ref_len = len(ref)
    alt_len = len(alt)
    snp_pos = rel_pos # relative position of the variant in the peak sequence (before adjusting for indels)
    input_len = peak_end - peak_start

    # get reference peak sequence first
    # peak_ref_seq = genome.fetch(chrom, peak_start, peak_end).upper()
    peak_ref_seq = cached_fetch(genome, chrom, peak_start, peak_end)

    # Check: construct REF peak sequence using variants info and check match with genome
    if check_ref_sequence:
        # extract REF variant sequence from genome and compare with input
        var_ref_seq = genome.fetch(chrom, var_pos - 1, var_pos - 1 + ref_len).upper()
        assert var_ref_seq == ref, print(f'var_ref_seq = {var_ref_seq}, ref = {ref}')

        # extract reference sequence from the peak sequence using relative coordinates
        var_ref_seq = peak_ref_seq[rel_pos:(rel_pos + ref_len)]
        assert var_ref_seq == ref, print(f'var_ref_seq = {var_ref_seq}, ref = {ref}')

        # get REF peak sequence by inserting reference allele (control matching ref sequence)
        peak_ref_seq_crl = peak_ref_seq[:rel_pos] + ref + peak_ref_seq[(rel_pos + ref_len):]
        assert peak_ref_seq_crl == peak_ref_seq, print(f'peak_ref_seq = {peak_ref_seq}, peak_ref_seq_crl = {peak_ref_seq_crl}')

    # find difference between ref/alt (indel)
    diff = ref_len - alt_len

    # handle deletions
    if diff > 0:
        # extend enough to allow space for all dels -> split deletion length to left/right
        left_extend = int(math.ceil(diff / 2))
        right_extend = int(math.floor(diff / 2))
        peak_start = peak_start - left_extend
        peak_end = peak_end + right_extend
        rel_pos = rel_pos + left_extend

        # get extended sequence with the alt-deleted allele
        peak_seq = cached_fetch(genome, chrom, peak_start, peak_end)
        peak_alt_seq = peak_seq[:rel_pos] + alt + peak_seq[(rel_pos + ref_len):]

        # check that alt peak sequence has correct length and sequence
        assert len(peak_alt_seq) == input_len, print('Length of alt sequence does not match input length!')
        assert peak_alt_seq[rel_pos:(rel_pos + alt_len)] == alt, print(f'peak_alt_seq = {peak_alt_seq}, alt = {alt}')

    # handle insertions
    elif diff < 0:
        # check rel_pos vs indel_len to determine whether and how much to remove from left/right
        diff = abs(diff)

        # if there isn't enough context left, trim from the other end
        if rel_pos - int(math.ceil(diff / 2)) < 25:
            peak_end = peak_end - diff
            # print(f'New peak_end = {peak_end}')
            # print('\n')
        elif input_len - rel_pos - int(math.ceil(diff / 2)) < 25:
            peak_start = peak_start + diff
            rel_pos = rel_pos - diff
            # print(f'New peak_start = {peak_start}')
            # print('\n')
        else:
            left_trim = int(math.ceil(diff / 2))
            right_trim = int(math.floor(diff / 2))
            peak_start = peak_start + left_trim
            peak_end = peak_end - right_trim
            rel_pos = rel_pos - left_trim

        # get trimmed sequence with the alt-inserted allele
        if peak_end - peak_start > 0:
            peak_seq = cached_fetch(genome, chrom, peak_start, peak_end)
            peak_alt_seq = peak_seq[:rel_pos] + alt + peak_seq[(rel_pos + ref_len):]
        else:
            peak_alt_seq = alt[:input_len]

        # check that alt peak sequence has correct length
        assert len(peak_alt_seq) == input_len, print('Length of alt sequence does not match input length!')
        assert peak_alt_seq[rel_pos:(rel_pos + alt_len)] == alt, print(f'peak_alt_seq = {peak_alt_seq}, alt = {alt}')

    # handle SNPs
    else:
        peak_seq = cached_fetch(genome, chrom, peak_start, peak_end)
        peak_alt_seq = peak_seq[:rel_pos] + alt + peak_seq[(rel_pos + 1):]

    return peak_ref_seq, peak_alt_seq, snp_pos

######

# register the genome
genome = crested.Genome(
    "/staging/leuven/res_00001/genomes/homo_sapiens/CHM13v2_maskedY_rCRS/fasta/chm13v2.0_maskedY_rCRS.fa", 
    "/staging/leuven/res_00001/genomes/homo_sapiens/CHM13v2_maskedY_rCRS/fasta/chm13v2.0_maskedY_rCRS.chrom.sizes"
)
crested.register_genome(
    genome
)

# CC deepPeak model
model_path = f"{model_dir}/deepPeak_CC_mean_finetuned/02.keras"
#model_path = "/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/deepPeak/deepPeak_CC_mean_finetuned/checkpoints/02.keras"
adata_path = f"{model_dir}/adata_cc_mean_norm.h5ad"

model_cc = keras.models.load_model(model_path)
adata_cc = ad.read_h5ad(adata_path)
adata_cc.obs_names 

model_path = "/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/peak_regression/deepPeak_SN_mean_finetuned/01.keras"
adata_path = "/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/peak_regression/adata_sn_mean_norm.h5ad"

adata_sn = ad.read_h5ad(adata_path)
datamodule_sn = crested.tl.data.AnnDataModule(
    adata_sn,
    batch_size=256,  # lower this if you encounter OOM errors
)
model_sn = keras.models.load_model(
    model_path,
    compile=False,  # Choose the basemodel with best validation loss/performance metrics
)

ct_sn = adata_sn.obs_names

models_dict = {"CC": model_cc, "SN": model_sn}
adata_dict = {"CC": adata_cc, "SN": adata_sn}

# get caQTL & ASCA variants with CREsted score > 0.5

data_dir = "/staging/leuven/stg_00090/ASA/analysis/Manuscript/Tables/"
df = pl.read_csv(f"{data_dir}/caQTL_ASCA_ISM_combined.csv.gz")

# select max scoring variant per peak, cell type and brain region, get high scoring CC microglia peaks
df_peaks = df.filter(
                pl.col("signif_caQTL_ASCA_no_effect_size_filter"),
                #pl.col("logfc").abs() > 0.,
            ).sort(
                pl.col("logfc").abs(), descending=True
            ).group_by(
                ["peak_id", "cell_type", "brain_region"]
            ).agg(pl.all().first()
            ).with_columns(
                (pl.when(pl.col("logfc") < 0).then(pl.lit("ref")).otherwise(pl.lit("alt"))).alias("high_allele")
            )

print(df_peaks.group_by("brain_region").agg(pl.col("variant_id").n_unique(), pl.col("peak_id").n_unique()))

# Contribution scores - max variant per peak

seqs_high_allele = {}
seqs_low_allele = {}

region_cts = df_peaks.select(
        "brain_region", "cell_type"
    ).unique(
    ).with_columns(
        pl.col("cell_type").replace(cell_types_corr).alias("cell_type_corr")
    )

for ct_row in region_cts.iter_rows():

    region = ct_row[0]
    cell_type = ct_row[1]
    cell_type_corr = ct_row[2] # cell type as in the model
    print(f"Processing {region} - {cell_type}")
    region_ct = f"{region};{cell_type_corr}"
    
    df_sub = df_peaks.filter(
        (pl.col("brain_region") == region) & (pl.col("cell_type") == cell_type)
    )

    # save variant list for the region/cell type
    var_dir = f"contribution_scores/{region}/variant_lists/"
    os.makedirs(var_dir, exist_ok=True)
    df_sub.write_csv(f"{var_dir}/{cell_type_corr}.csv")
    print(f"Number of variants to process: {df_sub.shape[0]}")

    seqs_high_allele[region_ct] = []
    seqs_low_allele[region_ct] = []

    # get sequences for all variants in the region/cell type
    for row in df_sub.iter_rows():

        var_id = row[3]
        peak_id = row[0]
        high_allele = row[len(row) - 1]
        seq_ref, seq_alt, var_pos = get_var_sequence(var_id, peak_id, genome)

        if high_allele == "ref":
            seqs_high_allele[region_ct].append(seq_ref)
            seqs_low_allele[region_ct].append(seq_alt)
        else:
            seqs_high_allele[region_ct].append(seq_alt)
            seqs_low_allele[region_ct].append(seq_ref)

    # get contribution scores for high and low alleles
    adata = adata_dict[region]
    model = models_dict[region]
    idx = adata.obs_names.get_indexer([cell_type_corr]).item() # index for the cell type in the model

    ctr_high = crested.tl.contribution_scores(input=seqs_high_allele[region_ct], 
                                target_idx = idx,
                                model=model, 
                                all_class_names = adata.obs_names.tolist(),
                                method='integrated_grad',
                                output_dir=f"contribution_scores/{region}/high_allele")
    
    ctr_low = crested.tl.contribution_scores(input=seqs_low_allele[region_ct], 
                            target_idx = idx,
                            model=model, 
                            all_class_names = adata.obs_names.tolist(),
                            method='integrated_grad',
                            output_dir=f"contribution_scores/{region}/low_allele")
