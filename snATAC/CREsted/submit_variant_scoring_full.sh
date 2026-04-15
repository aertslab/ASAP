#!/bin/bash

outdir=/staging/leuven/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/20250703_CREsted_variant_scoring # running for selected models


####################
#### CC models #####
####################

# variants in CC 
variants=/staging/leuven/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/20250127_CREsted_variant_scoring/variants_in_peaks.CC.txt.gz # all in peaks
peaks=/staging/leuven/stg_00090/ASA/analysis//2024_T2T_ATAC_analysis/CC/out/consensus_peak_calling/CC_CT2_consensus_regions.bed # set to "" for faster calculation

# Topic regression 

model_path="/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/topic_regression/cc_full_or20_20k_train_scaled/48.keras"
model_name="cc_full_deepTopicRegr_or20_20k"
sbatch sbatch_variant_scoring.sh $model_path $model_name $variants $outdir $peaks

# deepPeak model (mean signal, finetuned)

model_path="/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/deepPeak/deepPeak_CC_mean_finetuned/checkpoints/02.keras"
model_name="CC_deepPeak_mean_finetuned"
sbatch sbatch_variant_scoring.sh $model_path $model_name $variants $outdir $peaks


# DeepHumanBrain model

model_dir="/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/published_models"
model_name="DeepHumanBrain"
model_path="${model_dir}/${model_name}.keras"
sbatch sbatch_variant_scoring.sh "$model_path" "$model_name" "$variants" "$outdir" "$peaks"


####################
#### SN models #####
####################


# variants in SN
variants=/staging/leuven/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/20250127_CREsted_variant_scoring/variants_in_peaks.SN.txt.gz
peaks=/staging/leuven/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/SN/out/consensus_peak_calling/SN_CT2_consensus_regions.bed

# SN deepTopic regression model
model_path="/lustre1/project/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/SN/deeplearning/SN_topic_regression_model_t20_OR20/17.keras"
model_name="sn_full_deepTopicRegr_or20_20k"
sbatch sbatch_variant_scoring.sh $model_path $model_name $variants $outdir $peaks

# SN deepPeak model (mean signal, finetuned)
model_path="/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/deepPeak/deepPeak_SN_mean_finetuned/checkpoints/01.keras"
model_name="sn_deepPeak_mean_finetuned"
sbatch sbatch_variant_scoring.sh $model_path $model_name $variants $outdir $peaks

