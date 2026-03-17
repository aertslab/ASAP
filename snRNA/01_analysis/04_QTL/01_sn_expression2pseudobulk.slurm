#!/bin/bash
#SBATCH -N 1
#SBATCH -c 15
#SBATCH --mem 200G
#SBATCH --time 05:00:00
#SBATCH --cluster wice
#SBATCH -A XX
#SBATCH -p XX
#SBATCH --job-name="01_expression2bulk_sn"
#SBATCH --mail-type=ALL
#SBATCH --output=logs/R-%x.%j.out
#SBATCH --error=logs/R-%x.%j.err
#SBATCH --mail-user=XX

# activate env; change the path as appropriate
source ~/.bash_profile
conda_initialize /path/to/conda
conda activate /path/to/conda/envs/pertpy-env


python 01_sn_expression2pseudobulk.py