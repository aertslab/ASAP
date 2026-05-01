#!/bin/bash -l


module load HTSlib/1.20-GCC-10.3.0
module load SAMtools/1.20-GCC-10.3.0
module load pigz/2.8-GCCcore-10.3.0
module load mawk/1.3.4-20240622-GCCcore-10.3.0
module load ISA-L/2.31.0-GCCcore-10.3.0
module load zstd/1.5.6-GCCcore-10.3.0
module load zlib/1.3.1
module load zlib-ng/2.1.6-GCCcore-10.3.0

apptainer_run="apptainer run --cleanenv --no-home -B /lustre1,/staging,/data,${VSC_SCRATCH},/tmp"


bwa2_mem_threads=16
bwa2_mem_fasta="/staging/leuven/res_00001/genomes/homo_sapiens/CHM13v2_maskedY_rCRS/indexes/bwa_mem2/2.2.1/chm13v2.0_maskedY_rCRS.fa"



bwa_mem2_mapping () {

    local fastq_file="${1}"
    local sample_id="${2}"
    local output_dir="${3}"
    #local rg="${4}"

    local fastq_basename_prefix="${fastq_file##*/}"
    local fastq_basename_prefix="${fastq_basename_prefix%.fq.gz}"

    printf 'fastq_file="%s"\n' "${fastq_file}";
    printf 'output_dir="%s"\n\n' "${output_dir}";
    printf 'fastq_basename_prefix="%s"\n' "${fastq_basename_prefix}"

    mkdir -p "${output_dir}"
  
    # Get library ID from first FASTQ read: instrument ID, run number, flowcell ID and lane number.
    local library_id=$(gzip -d -c "${fastq_file}" | awk -F ':' -v OFS=':' '{ $1 = substr($1, 2); NF = 4; print $0; exit }')

    # Mapping.
    ${apptainer_run} \
        /staging/leuven/res_00001/software/vsn_containers/vibsinglecellnf-bwamaptools-bwa-mem2-2.2.1-zlibng-2.0.6_ghuls.img \
            bwa-mem2 mem \
                -p \
                -t ${bwa2_mem_threads} \
                -C \
                -R "@RG\tID:${fastq_basename_prefix}\tSM:${sample_id}\tLB:${library_id}__${sample_id}\tPL:ILLUMINA" \
                "${bwa2_mem_fasta}" \
                "${fastq_file}" \
           | samtools fixmate -u -m - - \
           | samtools view -@ 2 -b -o "${output_dir}/${fastq_basename_prefix}.bwa.out.fixmate.bam"
}



bwa_mem2_mapping "${@}"
