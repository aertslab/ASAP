#!/bin/bash -l

analysis_dir="/staging/leuven/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/20250107_allelic_analysis"

# Path to the Snakefile template
snakefile_template="Snakefile_by_cell_type"

# Base directory for donor-specific Snakemake runs
snakemake_dir="./snakemake_by_region_ct"

# full list of selected barcodes with cluster asignments
barcodes_path="/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/data/barcodes_per_donor/selected_barcodes_combined.tsv"

# extract regions and cell types
cut -f3,5 $barcodes_path | tail -n +2 |  sort -u > regions_cell_types.txt

# file with regions and cell types to exclude
exclude_file="exclude_list.txt"


# loop over regions and cell types by line

while IFS=$'\t' read -r region cell_type; do
  
  printf "Region: %-5s Cell Type: %s\n" "$region" "$cell_type"

  # Check if the region-cell_type pair exists in the exclude file
  if grep -q -w -P "$region\t$cell_type" "$exclude_file"; then
    echo "NOT submitted"
  else
    
    # Create the subdirectory for region_cell_type
    region_ct_dir="$snakemake_dir/$region/$cell_type"
    mkdir -p "$region_ct_dir"

    # Copy and rename the Snakefile template into the donor directory
    cp "$snakefile_template" "$region_ct_dir/Snakefile"

    # Create the config.yaml file with the required parameters
    config_file="$region_ct_dir/config.yaml"

    # Create the config.yaml file
    cat <<EOL > "$config_file"
brain_region: "$region"
cell_type: "$cell_type"
EOL

    #submit snakemake
    cp sbatch_snakemake.sh $region_ct_dir
    cd $region_ct_dir
    sbatch sbatch_snakemake.sh
    cd ../../../

  fi

done < "regions_cell_types.txt"


