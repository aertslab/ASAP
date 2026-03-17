#!/bin/bash
#SBATCH -N 1
#SBATCH -c 64
#SBATCH --mem 1500G
#SBATCH --time 02:00:00
#SBATCH --cluster wice
#SBATCH -A XX
#SBATCH -p XX
#SBATCH --job-name="results-all"
#SBATCH --mail-type=ALL
#SBATCH --output=logs/R-%x.%j.out
#SBATCH --error=logs/R-%x.%j.err
#SBATCH --mail-user=XX

source ~/.bash_profile
mambaforge_initialize /path/to/mambaforge

mamba activate polars_large
python 06_fdr_all.py
