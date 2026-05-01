#!/bin/bash 

#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem 80G
#SBATCH -A lp_cbd_stae
#SBATCH -p batch
#SBATCH --cluster wice 
#SBATCH --time 48:00:00
#SBATCH --output logs/wasp_%j_log
#SBATCH --error errors/wasp_%j_err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user="olga.sigalova@kuleuven.be"

source /lustre1/project/stg_00002/mambaforge/vsc35059/etc/profile.d/conda.sh
conda activate genetic_variation_0624
mkdir -p logs errors

snakemake --cores 16


