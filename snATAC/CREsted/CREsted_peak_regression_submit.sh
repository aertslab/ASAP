#!/bin/bash

# 1. CC model

# adata_path=/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/peak_regression/adata_cc_norm.h5ad - for predicting counts, in the loss function set multiplier=1
# model_name=deepPeak_CC_basemodel
# sbatch_script=/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/5_CREsted_peak_regression_sbatch.sh 

adata_path=/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/peak_regression/adata_cc_mean_norm.h5ad
model_name=deepPeak_CC_mean
sbatch_script=/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/5_CREsted_peak_regression_sbatch.sh 

sbatch $sbatch_script $adata_path $model_name

# 2. SN model
adata_path=/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/peak_regression/adata_sn_mean_norm.h5ad
model_name=deepPeak_SN_mean
sbatch_script=/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/5_CREsted_peak_regression_sbatch.sh

#sbatch $sbatch_script $adata_path $model_name
