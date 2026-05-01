#!/bin/bash -l

analysis_dir="/staging/leuven/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/20250107_allelic_analysis"

# Path to the Snakefile template
snakefile_template="Snakefile_template"

# Base directory for donor-specific Snakemake runs
snakemake_dir="./snakemake_by_donor"

# full list of selected barcodes with cluster asignments
barcodes_path="/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/data/barcodes_per_donor/selected_barcodes_combined.tsv"
donors_list=$(cut -f4 $barcodes_path | tail -n +2 |  sort -u)


# Loop over the donors

count=0

for donor_id in $donors_list; do

    # limit the number of donors processed
    ((count++))

    # Skip donors first N donors 
    # if ((count < 81)); then
    #     continue
    # fi

    # # Stop after N donors
    # if ((count > 100)); then
    #     break
    # fi

    
    echo "$donor_id"

    # Create the subdirectory for this donor
    donor_dir="$snakemake_dir/$donor_id"
    mkdir -p "$donor_dir"

    # Copy and rename the Snakefile template into the donor directory
    cp "$snakefile_template" "$donor_dir/Snakefile"

    # Create the config.yaml file with the required parameters
    config_file="$donor_dir/config.yaml"

    # Create the config.yaml file
    cat <<EOL > "$config_file"
analysis_dir: "$analysis_dir"
donor_id: "$donor_id"
EOL

    #submit snakemake
    cp sbatch_snakemake.sh $donor_dir
    cd $donor_dir
    sbatch sbatch_snakemake.sh
    cd ../../

done


