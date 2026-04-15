#!/bin/bash
#SBATCH -n 1
#SBATCH -c 24
#SBATCH --mem 256G
#SBATCH -A lp_big_wice_gpu
#SBATCH -p dedicated_big_gpu
#SBATCH --gpus-per-node 1
#SBATCH --cluster wice
#SBATCH --time 36:00:00
#SBATCH --output logs/deepPeak_%j_log
#SBATCH --error errors/deepPeak_%j_err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user="olga.sigalova@kuleuven.be"


# activate environment
source /staging/leuven/stg_00002/lcb/osiga/software/anaconda3/etc/profile.d/conda.sh
conda activate /lustre1/project/stg_00002/mambaforge/vsc35059/envs/crested_070525

# input data
adata_path=$1
model_name=$2
script=/staging/leuven/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/2_CREsted_models/5_CREsted_peak_regression_train.py

# run the script 
python $script $adata_path $model_name

