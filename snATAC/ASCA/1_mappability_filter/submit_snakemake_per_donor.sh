#!/bin/bash -l

wasp_dir="/staging/leuven/stg_00090/ASA/src/wasp_indel_v1.0/"
vcf_dir="/staging/leuven/stg_00090/ASA/analysis/WGS_results_cram"
bam_dir="/staging/leuven/stg_00090/ASA/analysis/20241010_per_patient_bams/per_patient_bams"
outdir="/staging/leuven/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/20241612_mappability_filter"

# Path to the Snakefile template
snakefile_template="Snakefile_template"

# Base directory for donor-specific Snakemake runs
snakemake_dir="./snakemake_by_donor"

# Brain regions
#brain_region=CC
brain_region=SN

# List of processed donors in CC
#processed_donors=('ASA_095' 'ASA_174' 'ASA_167' 'ASA_070' 'ASA_090' 'ASA_108' 'ASA_001' 'ASA_002' 'ASA_005')
processed_bams=$(ls /staging/leuven/stg_00090/ASA/analysis/20241612_mappability_filter/"$brain_region"/*sorted.full.bam)

# Initialize an empty array to store donor IDs
processed_donors=()

# Loop through the processed BAM files
for bam_file in $processed_bams; do
    # Extract the donor_id by removing directory and file extension
    donor_id=$(basename "$bam_file" | sed 's/\.mappability_filter\.sorted\.full\.bam//')
    processed_donors+=("$donor_id")
done

echo "Processed donors: ${processed_donors[@]}"


# Function to check if a donor_id is in the processed_donors list
is_processed() {
    local donor_id="$1"
    for processed in "${processed_donors[@]}"; do
        if [[ "$processed" == "$donor_id" ]]; then
            return 0 # Found
        fi
    done
    return 1 # Not found
}

# Loop through all BAM files in the directory

count=0

for bam_file in "$bam_dir"/"$brain_region"/*.bam; do

    # limit the number of donors processed
    ((count++))

    # Skip the loop iteration if count < 51 (skip jobs that are running already)
    # if ((count < 51)); then
    #     continue
    # fi
    # Submit for subset of donors
    # if ((count > 150)); then
    #     break
    # fi

    # Extract the donor ID (filename without extension)
    donor_id=$(basename "$bam_file" .bam)

    # Check if the donor_id is NOT in the processed_donors list
    if ! is_processed "$donor_id"; then
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
wasp_dir: "$wasp_dir"
bam_dir: "$bam_dir"
outdir: "$outdir"
vcf_dir: "$vcf_dir"
donor_id: "$donor_id"
brain_region: "$brain_region"
EOL

        #submit snakemake
        cp sbatch_snakemake.sh $donor_dir
        cp bwa_mem2_mapping.sh $donor_dir
        cd $donor_dir
        sbatch sbatch_snakemake.sh
        cd ../../

    fi
done

