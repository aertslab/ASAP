#!/bin/bash
#SBATCH -n 1
#SBATCH -c 20
#SBATCH --mem 700G 
#SBATCH -A lp_big_wice_gpu
#SBATCH -p dedicated_big_gpu
#SBATCH --gpus-per-node 1
#SBATCH --cluster wice
#SBATCH --time 12:00:00
#SBATCH --output logs/variants_%j_log
#SBATCH --error errors/variants_%j_err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user="olga.sigalova@kuleuven.be"

# run with 700G for DeepHumanBrain, other models are ok with 256G

# Prepare environment (install CREsted) 
source /staging/leuven/stg_00002/lcb/osiga/software/anaconda3/etc/profile.d/conda.sh
conda activate crested_env_070525



#### ---------------------------------------------------------------------------------------------------------------------------------------------------------- ####


## trained model path
model_path=$1

# path to save the results
model_name=$2

## variants and peaks to score
schema='peak' # either snp-centered or peak-based ('snp' or 'peak' or 'indel')
variants_path=$3  # format of [chr, start, end, pos, ref, alt, id]

## output directory
outdir=$4
output_path="${outdir}/${model_name}"
echo $model_name
echo "Output written to ${output_path}"


# path to the source code (variant_utils.py and variant.py are copied in this directory)
DEEP_PEAK_DIR=/staging/leuven/stg_00002/lcb/vkonst/Regulatory_Models/deep_peak_exp_crested/src

## reference genome path
ref_fasta_path='/staging/leuven/res_00001/genomes/homo_sapiens/CHM13v2_maskedY_rCRS/fasta/chm13v2.0_maskedY_rCRS.fa'
ref_chrom_path='/staging/leuven/res_00001/genomes/homo_sapiens/CHM13v2_maskedY_rCRS/fasta/chm13v2.0_maskedY_rCRS.chrom.sizes'

peaks_path=$5 # leave empty (peaks_path='') for quicker computing of shap values


# variant scoring options
batch_size=1024

num_shuffles=0 # number of permuted scores/user/leuven/350/vsc35059/ASA/analysis/analysis_Olga/3_T2T_analysis/3_indel_scoring/cc_full_deepPeak_balanced per SNP (not yet implemented)
num_peaks=400000 # number of peaks to use for peak percentile calculation ('all' or integer)

method='ism' # method to get variant scores ('shap' or 'ism')
normalize=false # normalize gradients (for shap)
attribution_length=10 # window around SNP to compute attribution scores (for shap)
forward_only=false # get forward or forward + reverse strand scores (for ism)

#### ---------------------------------------------------------------------------------------------------------------------------------------------------------- ####

# Start of pipeline
base_name=$(basename ${variants_path} .tsv)

# Check output paths
echo "Output path: $output_path"

# Create the appropriate directories
mkdir -p $output_path
cd $output_path

# Run the variant scoring

if [ "$forward_only" = true ]; then
    if [ "$normalize" = true ]; then
        variant_scoring_path='variant_scoring_'$method'_'$attribution_length'_fo_norm_'$base_name
        echo "Variant scoring path: $variant_scoring_path"
        mkdir -p $variant_scoring_path
        $DEEP_PEAK_DIR/variant.py \
         -v $variants_path \
         -sc $schema \
         -m $model_path \
         -g $ref_fasta_path \
         -chrs $ref_chrom_path \
         -o $variant_scoring_path \
         -p $peaks_path \
         -b $batch_size \
         -ns $num_shuffles \
         -np $num_peaks \
         -me $method \
         -al $attribution_length \
         -fo \
         -norm
    else
        variant_scoring_path='variant_scoring_'$method'_'$attribution_length'_fo_'$base_name
        echo "Variant scoring path: $variant_scoring_path"
        mkdir -p $variant_scoring_path
        $DEEP_PEAK_DIR/variant.py \
         -v $variants_path \
         -sc $schema \
         -m $model_path \
         -g $ref_fasta_path \
         -chrs $ref_chrom_path \
         -o $variant_scoring_path \
         -p $peaks_path \
         -b $batch_size \
         -ns $num_shuffles \
         -np $num_peaks \
         -me $method \
         -al $attribution_length \
         -fo
    fi
elif [ "$forward_only" = false ]; then
    if [ "$normalize" = true ]; then
        variant_scoring_path='variant_scoring_'$method'_'$attribution_length'_norm_'$base_name
        echo "Variant scoring path: $variant_scoring_path"
        mkdir -p $variant_scoring_path
        $DEEP_PEAK_DIR/variant.py \
         -v $variants_path \
         -sc $schema \
         -m $model_path \
         -g $ref_fasta_path \
         -chrs $ref_chrom_path \
         -o $variant_scoring_path \
         -p $peaks_path \
         -b $batch_size \
         -ns $num_shuffles \
         -np $num_peaks \
         -me $method \
         -al $attribution_length \
         -norm
    else
        variant_scoring_path='variant_scoring_'$method'_'$attribution_length'_'$base_name
        echo "Variant scoring path: $variant_scoring_path"
        mkdir -p $variant_scoring_path
        $DEEP_PEAK_DIR/variant.py \
         -v $variants_path \
         -sc $schema \
         -m $model_path \
         -g $ref_fasta_path \
         -chrs $ref_chrom_path \
         -o $variant_scoring_path \
         -p $peaks_path \
         -b $batch_size \
         -ns $num_shuffles \
         -np $num_peaks \
         -me $method \
         -al $attribution_length
    fi
else
    echo "Unrecognized parameter"
fi