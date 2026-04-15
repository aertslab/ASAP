module load mawk/1.3.4-20240123-GCCcore-10.3.0
module load BEDTools/2.31.1-GCC-10.3.0


#metadata_atac_vsn_preprocess_tsv='/lustre1/project/stg_00090/ASA/doc/metadata_atac_vsn_preprocess_20240126.tsv'
#metadata_atac_vsn_preprocess_tsv='metadata_atac_vsn_preprocess_25052023_with_updated_AT005.tsv'
#metadata_atac_vsn_preprocess_tsv='/lustre1/project/stg_00090/ASA/doc/metadata_atac_vsn_preprocess_20240527.tsv'
metadata_atac_vsn_preprocess_tsv='/lustre1/project/stg_00090/ASA/doc/metadata_atac_vsn_preprocess_20240624.tsv'

# All short BC IDs.
❯  hck -d '\t' -F short_bc_id "${metadata_atac_vsn_preprocess_tsv}" | sort | uniq | grep -v short_bc | tr '\n' ' '
AT001a AT001b AT001c AT002a AT002b AT003a AT003b AT003c AT003d AT004a AT004b AT004c AT004d AT005a AT005b AT005c AT005d AT006a AT006b AT006c AT006d AT007a AT007b AT007c AT007d AT009a AT009b AT009c AT009d AT010a AT010c AT011a AT011c AT012a AT012b AT013a AT013b AT013c AT014a AT014b AT014c AT014d AT014e AT016a AT016b AT016c AT016d AT017a AT017b AT018b AT018c AT019a AT019b AT019c AT019d AT019e AT019f AT020a AT020b AT020c AT020d AT021a AT021b AT021c AT021d AT021e AT021f AT021g AT021h AT022b AT023a AT023b AT023c AT023d AT023e AT023f AT023g AT023h AT024a AT025a AT025b AT025c AT025d AT025e AT025f AT025g AT025h AT026a AT026b AT026c AT026d AT026e AT026f AT026g AT026h MO001a MO001b MO001c MO002a MO002b MO002c MO002d MO003a MO003b MO004a MO004b MO007a MO007b MO008a MO010a MO010b MO011a MO011b MO012a MO012b MO012c MO012d MO013a MO014a MO014b MO014c MO014d MO015a MO015b MO015c MO015d MO016a MO016b MO016c MO016d MO017a MO017b MO017c MO017d MO017e MO017f MO018a MO018b MO018c MO018d MO018e MO018f MO018g MO018h MO019a MO019b MO019c MO019d MO020a MO020b MO020c MO020d MO021a MO021b MO021c MO021d MO021e MO021f MO021g MO021h MO022a MO022b MO022c

# List all short BC IDs which are not processed before or which have new sequencing runs.
❯  zet diff \
    <(hck -d '\t' -F short_bc_id "${metadata_atac_vsn_preprocess_tsv}" | sort | uniq -c | grep -v short_bc ) \
    <(hck -d '\t' -F short_bc_id /staging/leuven/stg_00090/ASA/doc/metadata_atac_vsn_preprocess_20240527.tsv | sort | uniq -c | grep -v short_bc)
     32 AT026a
     32 AT026b
     32 AT026c
     32 AT026d
     32 AT026e
     32 AT026f
     32 AT026g
     32 AT026h
      9 AT027a
      9 AT027b

❯  zet diff \
    <(hck -d '\t' -F short_bc_id "${metadata_atac_vsn_preprocess_tsv}" | sort | uniq -c | grep -v short_bc ) \
    <(hck -d '\t' -F short_bc_id /staging/leuven/stg_00090/ASA/doc/metadata_atac_vsn_preprocess_20240527.tsv | sort | uniq -c | grep -v short_bc) \
    | awk '{printf $2 " "}'
AT026a AT026b AT026c AT026d AT026e AT026f AT026g AT026h AT027a AT027b


short_bc_ids=( AT026a AT026b AT026c AT026d AT026e AT026f AT026g AT026h AT027a AT027b )


❯  hck -d '\t' -F technology -F techology_original "${metadata_atac_vsn_preprocess_tsv}" | tail -n +2 | sort -u
HYAv2   HyDrop-ATAC 96 ligation
MultiomeATAC    10x Multiome ATAC
ScaleBio        ScaleATAC
ScaleBioIH      ScaleATAC
ScaleBioIH      scaleATAC
ScaleBioIH-HYAv2        HyDrop-scaleATAC 96 ligation
ScaleBioIH6     ScaleATAC
scATAC  10x scATAC
scATAC  10x scATAC v1.1
scATAC  10x scATAC v2



short_bc_ids_all=( $(hck -d '\t' -F short_bc_id "${metadata_atac_vsn_preprocess_tsv}" | sort | uniq |grep -v short_bc | tr '\n' ' ') )


# Create some text files:
#   - Create file with commands for creating short BC based FASTQ symlinks.
#       - Create symlinks to all original FASTQ file, but rename them to have an unique name:
#           short barcode ID + unique "_PXXX" suffix + "_{R1,R2,BC}.fastq.gz"
#       - Get all necessary information from the metadata
#       - Natural sort on the following fields:
#         1. Short barcode ID
#         2. Seqeuncing run name (in YYYYMMDD_sequencer format)
#         3. FASTQ PE1 file name for PE1
#       - Create commands to create the correct symlinks.
#   - Create extended metadata ATAC TSV file.
hck \
        -d '\t' \
        -F short_bc_id \
        -F merged_sample_name \
        -F pool \
        -F technology \
        -F technology \
        -F sequencing_run_name \
        -F fastq_PE1_path \
        -F fastq_PE2_path \
        -F fastq_barcode_path \
        -F fastq_barcode_path \
        "${metadata_atac_vsn_preprocess_tsv}" \
    | mawk \
        -F '\t' \
        -v 'OFS=\t' \
        '
        {
            # Only process non-header lines.
            if (NR != 1) {
                # Reverse sequencing run name, so date comes before sequencer,
                # FASTQ files can be sorted by date.
                split($5, sequencing_run_name_array, "_");
                $5 = sequencing_run_name_array[2] "_" sequencing_run_name_array[1];
                print $0;
            }
        }' \
    | sort -t $'\t' -k 1,1V -k 5,5V -k 6,6V \
    | mawk \
        -F '\t' \
        -v 'OFS=\t' \
        -v 'create_short_bc_based_fastq_symlinks_filename=create_short_bc_based_fastq_symlinks.txt' \
        -v 'metadata_atac_vsn_preprocess_extended_filename=metadata_atac_vsn_preprocess_extended.tsv' \
        '
        {
            if (NR == 1) {
                # Write header for extended ATAC metadata TSV file.
                print "short_bc_id", "merged_id", "pool", "technology", "sequencing_run_name", "fastq_PE1", "fastq_PE2", "fastq_BC", "fastq_basename_prefix", "fastq_PE1_with_fastq_basename_prefix", "fastq_PE2_with_fastq_basename_prefix", "fastq_BC_with_fastq_basename_prefix" > metadata_atac_vsn_preprocess_extended_filename;
            }

            short_bc_id = $1;
            merged_sample_name = $2;
            pool = $3;
            technology = $4;
            sequencing_run = $5;
            fastq_PE1 = $6;
            fastq_PE2 = $7;
            fastq_BC = $8;

            short_bc_id_to_fastq_no[short_bc_id] += 1;

            if (  short_bc_id_to_fastq_no[short_bc_id] == 1 ) {
                # Create fastq/raw/${short_bc_id} directory for each short BC ID.
                print "\nmkdir -p fastq/raw/" short_bc_id > create_short_bc_based_fastq_symlinks_filename;
            }

            # Create FASTQ basename prefix: short barcode ID + unique "_PXXX" suffix + "_{R1,R2,BC}.fastq.gz"
            fastq_basename_prefix=sprintf("%s_P%03d", short_bc_id, short_bc_id_to_fastq_no[short_bc_id]);

            # Create symlinks to original FASTQ files, but rename them to:
            #   fastq/raw/${short_bc_id}/${short_bc_id} + unique "_PXXX" suffix + "_{R1,R2,BC}.fastq.gz"
            fastq_PE1_ln_cmd = sprintf("ln -s %s fastq/raw/%s/%s_R1.fastq.gz", fastq_PE1, short_bc_id, fastq_basename_prefix);
            fastq_PE2_ln_cmd = sprintf("ln -s %s fastq/raw/%s/%s_R2.fastq.gz", fastq_PE2, short_bc_id, fastq_basename_prefix);
            fastq_BC_ln_cmd = sprintf("ln -s %s fastq/raw/%s/%s_BC.fastq.gz", fastq_BC, short_bc_id, fastq_basename_prefix);

            # Create full paths to new FASTQ files for extended ATAC metadata TSV file.
            fastq_PE1_with_fastq_basename_prefix = sprintf("/staging/leuven/stg_00090/ASA/analysis/20240624_atac_preprocess/fastq/raw/%s/%s_R1.fastq.gz", short_bc_id, fastq_basename_prefix);
            fastq_PE2_with_fastq_basename_prefix = sprintf("/staging/leuven/stg_00090/ASA/analysis/20240624_atac_preprocess/fastq/raw/%s/%s_R2.fastq.gz", short_bc_id, fastq_basename_prefix);
            fastq_BC_with_fastq_basename_prefix = sprintf("/staging/leuven/stg_00090/ASA/analysis/20240624_atac_preprocess/fastq/raw/%s/%s_BC.fastq.gz", short_bc_id, fastq_basename_prefix);

            # Write symlink commands to symlinks txt file.
            print fastq_PE1_ln_cmd > create_short_bc_based_fastq_symlinks_filename;
            print fastq_PE2_ln_cmd > create_short_bc_based_fastq_symlinks_filename;
            print fastq_BC_ln_cmd > create_short_bc_based_fastq_symlinks_filename;

            # Write current entry to extended ATAC metadata TSV file.
            print $0, fastq_basename_prefix, fastq_PE1_with_fastq_basename_prefix, fastq_PE2_with_fastq_basename_prefix, fastq_BC_with_fastq_basename_prefix > metadata_atac_vsn_preprocess_extended_filename;
        }'



# Check create_short_bc_based_fastq_symlinks.txt and if OK create all FASTQ file symlinks:
bash create_short_bc_based_fastq_symlinks.txt


# Get raw all FASTQ read counts.
# /staging/leuven/stg_00002/lcb/ghuls/software/single_cell_toolkit/rust/target/x86_64-unknown-linux-gnu/release/fastq_count \
#     $(hck -d '\t' -F fastq_BC_with_fastq_basename_prefix metadata_atac_vsn_preprocess_extended.tsv | tail -n +2) \
#  > fastq_counts.tsv

# Get raw FASTQ read counts for new/extra FASTQ files.
/staging/leuven/stg_00002/lcb/ghuls/software/single_cell_toolkit/rust/target/x86_64-unknown-linux-gnu/release/fastq_count \
    $(hck -d '\t' -F fastq_BC_with_fastq_basename_prefix metadata_atac_vsn_preprocess_extended.tsv | tail -n +2 | rg $(echo ${short_bc_ids[@]} | tr ' ' '|')) \
  > fastq_counts.extra.tsv

# Concat previous all raw FASTQ read counts with extra FASTQ read counts.
cat /staging/leuven/stg_00090/ASA/analysis/20240527_atac_preprocess/fastq_counts.tsv fastq_counts.extra.tsv > fastq_counts.tsv



# Check 

# Get all short barcode IDs and their associated technology.
hck \
        -d '\t' \
        -F short_bc_id \
        -F technology \
        "${metadata_atac_vsn_preprocess_tsv}" \
  | tail -n +2 \
  | sort -u -t $'\t' -k 1,1V -k 2,2V


# Run barcode correction, trimming and mapping.
hck \
        -d '\t' \
        -F fastq_PE1_with_fastq_basename_prefix \
        -F fastq_PE2_with_fastq_basename_prefix \
        -F fastq_BC_with_fastq_basename_prefix \
        -F short_bc_id \
        -F technology \
        metadata_atac_vsn_preprocess_extended.tsv \
  | tail -n +2 \
  | mawk -F '\t' \
        -v output_dir="/staging/leuven/stg_00090/ASA/analysis/20240624_atac_preprocess" \
        '
        {
            printf "/staging/leuven/stg_00090/ASA/analysis/20240624_atac_preprocess/barcode_correction_and_trimming_and_mapping.sh %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, output_dir;
        }
        '


# Run barcode correction, trimming and mapping with hyperqueue.
hck \
        -d '\t' \
        -F fastq_basename_prefix \
        -F fastq_PE1_with_fastq_basename_prefix \
        -F fastq_PE2_with_fastq_basename_prefix \
        -F fastq_BC_with_fastq_basename_prefix \
        -F short_bc_id \
        -F technology \
        metadata_atac_vsn_preprocess_extended.tsv \
  | tail -n +2 \
  | mawk -F '\t' \
        -v output_dir="/staging/leuven/stg_00090/ASA/analysis/20240624_atac_preprocess" \
        '
        {
            printf "hq submit --cpus 16 --name bc_trim_map_%s --cwd hq_wd /staging/leuven/stg_00090/ASA/analysis/20240624_atac_preprocess/barcode_correction_and_trimming_and_mapping.sh %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, output_dir;
        }
        '


# Run BAM merging and fragments creation with hyperqueue.
hck \
        -d '\t' \
        -F merged_id \
        -F short_bc_id \
        metadata_atac_vsn_preprocess_extended.tsv \
  | tail -n +2 \
  | sort -k 1,1V -k 2,2V \
  | bedtools groupby -g 1 -c 2 -o distinct \
  | mawk -F '\t' \
        -v output_dir="/staging/leuven/stg_00090/ASA/analysis/20240410_atac_preprocess" \
        '
        {
            printf "hq submit --cpus 12 --name merge_and_frags_%s --cwd hq_wd /staging/leuven/stg_00090/ASA/analysis/20240624_atac_preprocess/merge_bams_and_create_fragments_files.sh %s %s %s\n", $1, $2, $1, output_dir;
        }
        '


# Run BAM merging for short_bc_id and fragments creation with hyperqueue.
hck \
        -d '\t' \
        -F short_bc_id \
        metadata_atac_vsn_preprocess_extended.tsv \
  | tail -n +2 \
  | sort -u -k 1,1V \
  | mawk -F '\t' \
        -v output_dir="/staging/leuven/stg_00090/ASA/analysis/20240624_atac_preprocess" \
        '
        {
            printf "hq submit --cpus 12 --name merge_and_frags_%s --cwd hq_wd /staging/leuven/stg_00090/ASA/analysis/20240624_atac_preprocess/merge_bams_per_short_bc_id_and_create_fragments_files.sh %s %s\n", $1, $1, output_dir;
        }
        ' > merge_bams_per_short_bc_id_and_create_fragments_files.queue.txt



# Make hardlinks to previous run correct samples:
find /staging/leuven/stg_00090/ASA/analysis/20240126_atac_preprocess/fastq/barcode_correction_and_trimming/ \
    -type f \
    -name  '*_R*.fastq.gz*' \
    -o -name '*_R*_val_*.fq.gz*' \
    -o -name '*_corrected_bc.*' \
  | sort -V \
  | mawk '
    {
        old_filename = $1;
        new_filename = old_filename;

        # Fix new filename.
        sub("20240126_atac_preprocess", "20240410_atac_preprocess", new_filename);

        print "if [ ! -e \"" new_filename "\" ] ; then";
        print "    mkdir -p $(dirname \"" new_filename "\");";
        print "    ln \"" old_filename "\" \"" new_filename "\";";
        print "fi\n";
    }' > create_barcode_correction_and_trimming_hardlinks.txt

find /staging/leuven/stg_00090/ASA/analysis/20240126_atac_preprocess/mapping/per_short_bc_id/ \
    -type f \
    -name '*.bwa.out.fixmate.possorted.bam*' \
  | mawk '
    {
        old_filename = $1;
        new_filename = old_filename;

        # Fix new filename.
        sub("20240126_atac_preprocess", "20240410_atac_preprocess", new_filename);

        print "if [ ! -e \"" new_filename "\" ] ; then";
        print "    mkdir -p $(dirname \"" new_filename "\");";
        print "    ln \"" old_filename "\" \"" new_filename "\";";
        print "fi\n";
    }' > create_per_short_bc_id_hardlinks.txt

find /staging/leuven/stg_00090/ASA/analysis/20230110_atac_preprocess/mapping/per_short_bc_id_merged/ \
    -type f \
    -name '*.bwa.out.fixmate.possorted.bam*' \
    -o -name '*.fragments.raw.tsv.gz*' \
    -o -name '*.barcard.overlap.tsv' \
  | mawk '
    {
        old_filename = $1;
        new_filename = old_filename;

        # Fix new filename.
        sub("20230110_atac_preprocess", "20230525_atac_preprocess", new_filename);

        print "if [ ! -e \"" new_filename "\" ] ; then";
        print "    mkdir -p $(dirname \"" new_filename "\");";
        print "    ln \"" old_filename "\" \"" new_filename "\";";
        print "fi\n";
    }' > create_per_short_bc_id_merged_hardlinks.txt


# Run BAM merging for short_bc_id and fragments creation with hyperqueue.
hck \
        -d '\t' \
        -F short_bc_id \
        /staging/leuven/stg_00090/ASA/analysis/20230525_atac_preprocess/metadata_atac_vsn_preprocess_extended.tsv \
  | tail -n +2 \
  | sort -u -k 1,1V \
  | mawk -F '\t' \
        -v output_dir="/staging/leuven/stg_00090/ASA/analysis/20230525_atac_preprocess" \
        '
        {
            printf "hq submit --cpus 6 --name mapping_stats_%s /staging/leuven/stg_00090/ASA/analysis/20230525_atac_preprocess/mapping_stats_per_short_bc_id.sh %s %s\n", $1, $1, output_dir;
        }
        ' > mapping_stats_per_short_bc_id_merged.queue.txt
