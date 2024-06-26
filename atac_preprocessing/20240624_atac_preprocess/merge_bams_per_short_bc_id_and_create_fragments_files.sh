#!/bin/bash -l


module load HTSlib/1.20-GCC-10.3.0
module load SAMtools/1.20-GCC-10.3.0
module load pigz/2.7-GCCcore-10.3.0
module load mawk/1.3.4-20240123-GCCcore-10.3.0
module load ISA-L/2.30.0-GCCcore-10.3.0
module load zstd/1.5.5-GCCcore-10.3.0
module load zlib/1.2.13
module load zlib-ng/2.1.6-GCCcore-10.3.0


apptainer_run="apptainer run --cleanenv --no-home -B /lustre1,/staging,/data,${VSC_SCRATCH},/tmp"


trim_galore_threads=8
bwa2_mem_threads=16
bwa2_mem_fasta="/staging/leuven/res_00001/genomes/homo_sapiens/CHM13v2_maskedY_rCRS/indexes/bwa_mem2/2.2.1/chm13v2.0_maskedY_rCRS.fa"



merge_bams () {
    local short_bc_id="${1}";
    local output_dir="${2}";

    local mapping_per_short_bc_id_dir="${output_dir}/mapping/per_short_bc_id";
    local mapping_per_short_bc_id_merged_dir="${output_dir}/mapping/per_short_bc_id_merged/${short_bc_id}";

    if [ "${short_bc_id/,/}" = "${short_bc_id}" ] ; then
        short_bc_ids_bam_files=$(echo "${mapping_per_short_bc_id_dir}/"${short_bc_id}/*.bwa.out.fixmate.possorted.bam);
    else
        short_bc_ids_bam_files=$(eval echo $(echo "${mapping_per_short_bc_id_dir}/"{${short_bc_id}})/*.bwa.out.fixmate.possorted.bam);
    fi

    mkdir -p "${mapping_per_short_bc_id_merged_dir}"

    samtools merge \
        -@ 8 \
        --write-index \
        -o "${mapping_per_short_bc_id_merged_dir}/${short_bc_id}.bwa.out.fixmate.possorted.bam##idx##${mapping_per_short_bc_id_merged_dir}/${short_bc_id}.bwa.out.fixmate.possorted.bam.bai" \
        ${short_bc_ids_bam_files};
}



mapping_stats () {
    local short_bc_id="${1}";
    local output_dir="${2}";

    local mapping_per_short_bc_id_merged_dir="${output_dir}/mapping/per_short_bc_id_merged/${short_bc_id}";
    local mapping_stats_per_short_bc_id_merged_dir="${output_dir}/mapping_stats/per_short_bc_id_merged";

    mkdir -p "${mapping_stats_per_short_bc_id_merged_dir}"

    /staging/leuven/stg_00002/lcb/ghuls/software/single_cell_toolkit/mapping_stats.sh \
        "${short_bc_id}" \
        "${mapping_per_short_bc_id_merged_dir}/${short_bc_id}.bwa.out.fixmate.possorted.bam" \
        "${mapping_stats_per_short_bc_id_merged_dir}";
}



create_fragments_file () {
    local short_bc_id="${1}";
    local output_dir="${2}";

    local singlecelltoolkit_apptainer_run="${apptainer_run} /staging/leuven/res_00001/software/vsn_containers/vibsinglecellnf-singlecelltoolkit-2024-04-09-62429e9.img";

    local mapping_per_short_bc_id_merged_dir="${output_dir}/mapping/per_short_bc_id_merged/${short_bc_id}";

    #${singlecelltoolkit_apptainer_run} \
    #    create_fragments_file \
    ${singlecelltoolkit_apptainer_run} \
        create_fragments_file \
            --bam "${mapping_per_short_bc_id_merged_dir}/${short_bc_id}.bwa.out.fixmate.possorted.bam" \
            --fragments "${mapping_per_short_bc_id_merged_dir}/${short_bc_id}.fragments.raw.tsv.gz";

    tabix -p bed "${mapping_per_short_bc_id_merged_dir}/${short_bc_id}.fragments.raw.tsv.gz";

    local chromosome_regex='^(chr)?([0-9]+|[XY])\$';

    ${singlecelltoolkit_apptainer_run} \
        calculate_jaccard_index_cbs.py \
            -i "${mapping_per_short_bc_id_merged_dir}/${short_bc_id}.fragments.raw.tsv.gz" \
            -o "${mapping_per_short_bc_id_merged_dir}/${short_bc_id}.barcard.overlap.tsv" \
            -t 1000 \
            -c "${chromosome_regex}";
}



merge_bams_and_create_fragments_files () {
    local short_bc_id="${1}";
    local output_dir="${2}";

    merge_bams "${short_bc_id}" "${output_dir}";
    mapping_stats "${short_bc_id}" "${output_dir}";

    create_fragments_file "${short_bc_id}" "${output_dir}";
}



merge_bams_and_create_fragments_files "${@}";
