#!/bin/bash -l


module load HTSlib/1.19.1-GCC-10.3.0
module load SAMtools/1.19.1-GCC-10.3.0
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


trim_galore () {
    local fastq_PE1_file="${1}";
    local fastq_PE2_file="${2}";
    local output_dir="${3}";

    local fastq_PE1_trimmed_file="${output_dir}/$(basename "${fastq_PE1_file%.fq.gz}_val_1.fq.gz")";
    local fastq_PE2_trimmed_out="${output_dir}/$(basename "${fastq_PE2_file%.fq.gz}_val_2.fq.gz")";

    ${apptainer_run} \
        /staging/leuven/res_00001/software/vsn_containers/vibsinglecellnf-trimgalore-trimgalore-0.6.7-cutadapt-4.2.img \
            trim_galore \
                -j "${trim_galore_threads}" \
                --paired \
                --gzip \
                -o "${output_dir}" \
                "${fastq_PE1_file}" \
                "${fastq_PE2_file}"
}



correct_10x_atac_barcode_from_fastq () {
    local fastq_BC_file="${1}";
    local corrected_BC_file="${2}";
    local bc_suffix="${3}";

    /staging/leuven/stg_00002/lcb/ghuls/software/single_cell_toolkit/correct_barcode_from_fastq.sh \
        /staging/leuven/res_00001/barcodes/cellranger_atac.737K-cratac-v1.txt.gz \
        "false" \
        "${fastq_BC_file}" \
        "${corrected_BC_file}" \
        "${bc_suffix}"
}



correct_hydrop_atac_barcode_from_fastq () {
    local fastq_BC_file="${1}";
    local corrected_BC_file="${2}";
    local bc_suffix="${3}";

    /staging/leuven/stg_00002/lcb/ghuls/software/single_cell_toolkit/extract_and_correct_hydrop_atac_barcode_from_fastq.sh \
        /staging/leuven/res_00001/barcodes/HyDrop_v2.txt \
        "${fastq_BC_file}" \
        "${corrected_BC_file}" \
        "${bc_suffix}"
}


correct_10x_multiome_atac_barcode_from_fastq () {
    local fastq_BC_file="${1}";
    local corrected_BC_file="${2}";
    local bc_suffix="${3}";

    /staging/leuven/stg_00002/lcb/ghuls/software/single_cell_toolkit/correct_barcode_from_fastq.sh \
        /staging/leuven/res_00001/barcodes/cellranger_arc_atac.737K-arc-v1.txt.gz \
        /staging/leuven/res_00001/barcodes/cellranger_arc_rna.737K-arc-v1.txt.gz \
        "${fastq_BC_file}" \
        "${corrected_BC_file}" \
        "${bc_suffix}"
}


extract_and_correct_scalebio_atac_barcode_from_fastq () {
    local fastq_BC_file="${1}";
    local corrected_BC_file="${2}";
    local scalebio_bc_type="${3}";
    local bc_suffix="${4}";

    case "${scalebio_bc_type}" in
        scalebio|scalebioih|scalebioih6)
            local tenx_or_hydrop_atac_bc_whitelist_file='/staging/leuven/res_00001/barcodes/cellranger_atac.737K-cratac-v1.txt.gz';;
        scalebioih-hyav2)
            local tenx_or_hydrop_atac_bc_whitelist_file='/staging/leuven/res_00001/barcodes/HyDrop_v2.txt';
            scalebio_bc_type='scalebioih';;
        *)
            printf 'Error: Unsupported scalebio BC type "%s".\n' "${scalebio_bc_type}";
            return 1;;
    esac

    /staging/leuven/stg_00002/lcb/ghuls/software/single_cell_toolkit/extract_and_correct_scalebio_atac_barcode_from_fastq.sh \
        "${tenx_or_hydrop_atac_bc_whitelist_file}" \
        "${fastq_BC_file}" \
        "${corrected_BC_file}" \
        "${scalebio_bc_type}" \
        false \
        "${bc_suffix}"
}


barcode_correction_and_trimming_and_mapping () {
    local fastq_PE1_file="${1}";
    local fastq_PE2_file="${2}";
    local fastq_BC_file="${3}";
    local short_bc_id="${4}";
    local technology="${5}";
    local output_dir="${6}";

    if [ ${#@} -ne 6 ] ; then
        printf 'Usage: barcode_correction_and_trimming_and_mapping fastq_PE1_file fastq_PE2_file fastq_BC_file short_bc_id technology output_dir\n';
        return 1;
    fi


    local barcode_correction_and_trimming_short_bc_id_dir="${output_dir}/fastq/barcode_correction_and_trimming/${short_bc_id}"
    local mapping_per_short_bc_id_dir="${output_dir}/mapping/per_short_bc_id/${short_bc_id}"

    local fastq_basename_prefix="${fastq_PE1_file##*/}"
    local fastq_basename_prefix="${fastq_basename_prefix%_R1.fastq.gz}"

    local corrected_BC_file="${output_dir}/fastq/barcode_correction_and_trimming/${short_bc_id}/${fastq_basename_prefix}_corrected_bc.zst";

    local fastq_PE1_trim_galore_file="${output_dir}/fastq/barcode_correction_and_trimming/${short_bc_id}/${fastq_basename_prefix}_R1_val_1.fq.gz";
    local fastq_PE2_trim_galore_file="${output_dir}/fastq/barcode_correction_and_trimming/${short_bc_id}/${fastq_basename_prefix}_R2_val_2.fq.gz";

    local fastq_PE1_with_corrected_BC_file="${output_dir}/fastq/barcode_correction_and_trimming/${short_bc_id}/${fastq_basename_prefix}_R1.fastq.gz";
    local fastq_PE2_with_corrected_BC_file="${output_dir}/fastq/barcode_correction_and_trimming/${short_bc_id}/${fastq_basename_prefix}_R2.fastq.gz";


    printf 'fastq_PE1_file="%s"\n' "${fastq_PE1_file}";
    printf 'fastq_PE2_file="%s"\n' "${fastq_PE2_file}";
    printf 'fastq_BC_file="%s"\n' "${fastq_BC_file}";
    printf 'short_bc_id="%s"\n' "${short_bc_id}";
    printf 'technology="%s"\n' "${technology}";
    printf 'output_dir="%s"\n\n' "${output_dir}";

    printf 'barcode_correction_and_trimming_short_bc_id_dir="%s"\n' "${barcode_correction_and_trimming_short_bc_id_dir}";
    printf 'mapping_per_short_bc_id_dir="%s"\n\n' "${mapping_per_short_bc_id_dir}";

    printf 'fastq_basename_prefix="%s"\n' "${fastq_basename_prefix}";
    printf 'corrected_BC_file="%s"\n' "${corrected_BC_file}";
    printf 'fastq_PE1_trim_galore_file="%s"\n' "${fastq_PE1_trim_galore_file}";
    printf 'fastq_PE2_trim_galore_file="%s"\n' "${fastq_PE2_trim_galore_file}";
    printf 'fastq_PE1_with_corrected_BC_file="%s"\n' "${fastq_PE1_with_corrected_BC_file}";
    printf 'fastq_PE2_with_corrected_BC_file="%s"\n\n' "${fastq_PE2_with_corrected_BC_file}";

    # return 0;

#    if [ -e "${fastq_PE2_with_corrected_BC_file}" ] ; then
#        printf 'Warning: FASTQ PE2 file "%s" exist. Assume "%s" was processed correctly.\n' \
#            "${fastq_PE2_with_corrected_BC_file}" \
#            "${fastq_basename_prefix}";
#        return 0;
#    fi


    if samtools quickcheck "${mapping_per_short_bc_id_dir}/${fastq_basename_prefix}.bwa.out.fixmate.possorted.bam" ; then
    #if [ -e "${mapping_per_short_bc_id_dir}/${fastq_basename_prefix}.bwa.out.fixmate.possorted.bam.bai" ] ; then
        printf 'Warning: BAM file "%s" is complete. Assume "%s" was processed correctly.\n' \
            "${mapping_per_short_bc_id_dir}/${fastq_basename_prefix}.bwa.out.fixmate.possorted.bam" \
            "${fastq_basename_prefix}";
        return 0;
    fi

    mkdir -p "${barcode_correction_and_trimming_short_bc_id_dir}" "${mapping_per_short_bc_id_dir}";

    # Check technology (after converting to lowercase) and run correct barcode correction script.
    case "${technology,,}" in
        multiomeatac)
            correct_10x_multiome_atac_barcode_from_fastq \
                "${fastq_BC_file}" \
                "${corrected_BC_file}" \
                "${short_bc_id}";;
        scalebio|scalebioih|scalebioih6|scalebioih-hyav2)
            extract_and_correct_scalebio_atac_barcode_from_fastq \
                "${fastq_BC_file}" \
                "${corrected_BC_file}" \
                "${technology,,}" \
                "${short_bc_id}";;
        hyav2)
            correct_hydrop_atac_barcode_from_fastq \
                "${fastq_BC_file}" \
                "${corrected_BC_file}" \
                "${short_bc_id}";;
        scatac)
            correct_10x_atac_barcode_from_fastq \
                "${fastq_BC_file}" \
                "${corrected_BC_file}" \
                "${short_bc_id}";;
        *)
            printf 'Error: Unsupported technology "%s".\n' "${technology}";
            return 1;;
    esac

    if [ ! -e "${corrected_BC_file}" ] ; then
        printf 'Error: Could not create corrected BC filename: "%s"\n' "${corrected_BC_file}";
        return 1;
    fi

    #if [ ! -e "${fastq_PE2_with_corrected_BC_file}" ] ; then

    ${apptainer_run} \
        /staging/leuven/res_00001/software/vsn_containers/vibsinglecellnf-trimgalore-trimgalore-0.6.7-cutadapt-4.2.img \
            trim_galore \
                -j "${trim_galore_threads}" \
                --paired \
                --gzip \
                -o "${barcode_correction_and_trimming_short_bc_id_dir}" \
                "${fastq_PE1_file}" \
                "${fastq_PE2_file}"

    /staging/leuven/stg_00002/lcb/ghuls/software/single_cell_toolkit/add_corrected_barcode_to_read_name.sh \
        "${fastq_PE1_trim_galore_file}" \
        "${fastq_PE2_trim_galore_file}" \
        "${corrected_BC_file}" \
        "${barcode_correction_and_trimming_short_bc_id_dir}/${fastq_basename_prefix}"

    #fi

    # Get library ID from first FASTQ read: instrument ID, run number, flowcell ID and lane number.
    local library_id=$(gzip -d -c "${fastq_PE1_with_corrected_BC_file}" | awk -F ':' -v OFS=':' '{ $1 = substr($1, 2); NF = 4; print $0; exit }');

    # Mapping.
    ${apptainer_run} \
        /staging/leuven/res_00001/software/vsn_containers/vibsinglecellnf-bwamaptools-bwa-mem2-2.2.1-zlibng-2.0.6_ghuls.img \
            bwa-mem2 mem \
                -t ${bwa2_mem_threads} \
                -C \
                -R "@RG\tID:${fastq_basename_prefix}\tSM:${short_bc_id}\tLB:${library_id}__${short_bc_id}\tPL:ILLUMINA" \
                "${bwa2_mem_fasta}" \
                "${fastq_PE1_with_corrected_BC_file}" \
                "${fastq_PE2_with_corrected_BC_file}" \
            | samtools fixmate -u -m -O bam - - \
            | samtools sort \
                -@ 2 \
                -m 2G \
                -O bam \
                --write-index \
                -T "/tmp/${fastq_basename_prefix}.bwa.out.fixmate.possorted.TMP" \
                -o "${mapping_per_short_bc_id_dir}/${fastq_basename_prefix}.bwa.out.fixmate.possorted.bam##idx##${mapping_per_short_bc_id_dir}/${fastq_basename_prefix}.bwa.out.fixmate.possorted.bam.bai" \
                -
}



barcode_correction_and_trimming_and_mapping "${@}"
