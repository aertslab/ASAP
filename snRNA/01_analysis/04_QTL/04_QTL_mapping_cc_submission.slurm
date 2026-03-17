#!/bin/bash
#SBATCH -N 1
#SBATCH -c 15
#SBATCH --mem 100G
#SBATCH --time 20:00:00
#SBATCH --cluster wice
#SBATCH -A XX_gpu
#SBATCH -p XX_gpu
#SBATCH --gpus-per-node 1
#SBATCH --job-name="04_tesnorQTL_cc"
#SBATCH --mail-type=ALL
#SBATCH --output=logs/R-%x.%j.out
#SBATCH --error=logs/R-%x.%j.err
#SBATCH --mail-user=XX
#SBATCH --array=1-2

names=("04_QTL_mapping_cc.py" "04_QTL_mapping_indel_cc.py")
name=${names[${SLURM_ARRAY_TASK_ID}]}

singularity exec --nv \
  -B /staging,/lustre1,/data,/vsc-hard-mounts,/local_scratch,/scratch,/user,/apps \
  /software/tensorqtl_1.0.10-v3.sif \
  python3 ${name}
