#! bin/bash

source /lustre1/project/stg_00002/mambaforge/vsc35059/etc/profile.d/conda.sh
conda activate genetic_variation_0624

snakemake --cores 4 --rerun-incomplete