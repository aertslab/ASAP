import polars as pl
import pandas as pd
import glob
import re
import os
import yaml
import pysam

# Define analysis directory and BCF path
analysis_dir = "/staging/leuven/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/20250505_allelic_analysis"
bcf_path = f"{analysis_dir}/WGS_chm13_BCFs.missing_to_ref.norm.subset_peaks_exended2500bp.vcf.gz"
bcf = pysam.VariantFile(bcf_path)

# Get contig names and samples from the BCF
chroms = [x for x in bcf.header.contigs]
samples = [x for x in bcf.header.samples]

# Make output directory
outdir = os.path.join(analysis_dir, "phasing_blocks")
os.makedirs(f"{outdir}/per_chrom", exist_ok=True)

# Parse BCF file and extract phasing information
for chrom in chroms:
    print(chrom)
    phasing_dict = {}
    i = 0

    for record in bcf.fetch(contig=chrom):
        for sample in record.samples:
            try:
                phasing_block = record.samples[sample]['PS']
                if phasing_block:
                    entry = [record.chrom, record.pos, sample, phasing_block]
                    phasing_dict.setdefault(sample, []).append(entry)
            except:
                continue

        if i % 200000 == 0:
            print(i)
        i += 1

    # Write phasing blocks to file
    for sample in samples:
        try:
            df = pd.DataFrame(phasing_dict[sample], columns=['chrom', 'pos', 'sample', 'phasing_block'])
            df = pl.from_pandas(df)
            (
                df
                .select("chrom", "pos", "phasing_block")
                .group_by(["chrom", "phasing_block"])
                .agg(pl.col("pos").min().alias("start"), pl.col("pos").max().alias("end"))
                .select(["chrom", "start", "end", "phasing_block"])
                .sort("chrom", "start")
                .write_csv(f"{outdir}/per_chrom/pb_{sample}_{chrom}.bed", separator="\t", include_header=False)
            )
        except:
            continue

# Extract sample names from per-chromosome BEDs
sample_names = []
for file in glob.glob(f"{outdir}/per_chrom/*_chr*.bed"):
    match = re.search(r"pb_(ASA_[0-9]+)_chr*", file)
    if match:
        sample_names.append(match.group(1))

phased_samples = set(sample_names)

# Merge per-chromosome BED files for each sample
os.makedirs(f"{outdir}/merged", exist_ok=True)
for sample in phased_samples:
    chrom_files = glob.glob(f"{outdir}/per_chrom/pb_{sample}_chr*.bed")
    with open(f"{outdir}/merged/pb_{sample}.bed", "w") as outfile:
        for fname in chrom_files:
            with open(fname) as infile:
                outfile.write(infile.read())