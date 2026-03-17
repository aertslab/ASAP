#!/bin/bash
#SBATCH -N 1
#SBATCH -c 15
#SBATCH --mem 100G
#SBATCH --time 20:00:00
#SBATCH --cluster wice
#SBATCH -A XX
#SBATCH -p XX
#SBATCH --gpus-per-node 1
#SBATCH --job-name="02_tesnorQTL_sn_SNP_indel_sv"
#SBATCH --mail-type=ALL
#SBATCH --output=logs/R-%x.%j.out
#SBATCH --error=logs/R-%x.%j.err
#SBATCH --mail-user=XX

# change paths to the script and containter as necessary
singularity exec --nv \
  -B /staging,/lustre1,/data,/vsc-hard-mounts,/local_scratch,/scratch,/user,/apps \
  /software/tensorqtl_1.0.10-v3.sif \
  python3 02_QTL_mapping_SNP_sn.py
