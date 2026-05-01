#!/bin/bash 

#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem 80G
#SBATCH -A lp_cbd_stae
#SBATCH -p batch
#SBATCH --cluster wice 
#SBATCH --time 3:00:00
#SBATCH --output logs/as_counts_%j_log
#SBATCH --error errors/as_counts_%j_err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user="olga.sigalova@kuleuven.be"

source /lustre1/project/stg_00002/mambaforge/vsc35059/etc/profile.d/conda.sh
conda activate genetic_variation_0624
mkdir -p logs errors

# only use this if sure that all files generated earlier are correct (to execute only new rules in the snakemake run below)
#snakemake --touch --cores 16

# run the full pipeline (or the new parts)
snakemake --cores 16

# rerun specific rule
# snakemake --cores 16 --forcerun get_target_region_counts 