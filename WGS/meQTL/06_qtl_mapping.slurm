#!/bin/bash
#SBATCH -N 1
#SBATCH -c 20
#SBATCH --mem 250G
#SBATCH --time 20:00:00
#SBATCH --cluster wice
#SBATCH -A XX_gpu_h100
#SBATCH -p XX_gpu_h100
#SBATCH --gpus-per-node 1
#SBATCH --job-name="06_qtl_mapping"
#SBATCH --mail-type=ALL
#SBATCH --output=logs/R-%x.%j.out
#SBATCH --error=logs/R-%x.%j.err
#SBATCH --mail-user=XX

# upade the path to the tensorqtl singularity container and the script as needed 
singularity exec --nv \
  -B /staging,/lustre1,/data,/vsc-hard-mounts,/local_scratch,/scratch,/user,/apps \
  /software/tensorqtl_1.0.10-v2.sif \
  python3 06_qtl_mapping.py
