#!/bin/bash
#SBATCH -n 1
#SBATCH -c 36
#SBATCH --time=10:00:00
#SBATCH -A XX
#SBATCH -p XX
#SBATCH --cluster wice
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XX
#SBATCH --output=logs/R-%x.%j.out
#SBATCH --error=logs/R-%x.%j.err
#SBATCH --job-name=sniffles_2.6_pop_call
# script to combine sniffles snf files to single multisample vcf

# change paths as needed
SINGULARITY_SNIFFLES="/software/singularity/sniffles_2.6.0.sif" 
# tsv file with snf files to do population calling on
# to get TSV with IDs
# awk -F'/' '{print $0 "\t" $(NF-2)}' ASAP_cleaned.tsv > ASAP_cleaned_IDs.tsv 
INPUT="ASAP_cleaned_IDs.tsv"
current_date=$(date "+%Y%m%d")

# run sniffles
singularity run -B /staging/leuven/ \
$SINGULARITY_SNIFFLES sniffles --input $INPUT \
--combine-max-inmemory-results 200 --combine-pctseq 0.5  --threads 36 --vcf ${current_date}_sniffles_2.6_ASAP.pct0.5.vcf.gz
