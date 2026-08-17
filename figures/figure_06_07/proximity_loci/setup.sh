#!/usr/bin/env bash
set -euo pipefail

DONOR_BCF=""
PASSTHROUGH_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --donor-bcf)
            DONOR_BCF="$2"
            shift 2
            ;;
        *)
            PASSTHROUGH_ARGS+=("$1")
            shift
            ;;
    esac
done

echo
echo " Running proximity (LD) block setup."
echo " This can take half an hour if all the steps have to be executed."
echo " Run in a compute session if available."
echo " Downloads plink2, lead variants and summary from GP2 GWAS,"
echo " 1KG T2T bcf, runs pgen creation and r2 calc using plink2"
echo

ENV_NAME="PROX"
ENV_FILE="prox.yml"

# Create the env if it doesn't already exist
echo " [setup] Checking PROX env..."
if ! conda env list | grep -qE "^\s*${ENV_NAME}\s"; then
    printf "\n [env] Creating conda env '${ENV_NAME}' from ${ENV_FILE}\n"
    conda env create -f "${ENV_FILE}"
fi

# Make the folder structure
echo " [setup] Creating folders..."
mkdir -p 1_plink
mkdir -p 2_source
mkdir -p 3_variants
mkdir -p 4_regions
mkdir -p 5_figs
mkdir -p 6_enrichment

# Download T2T fasta
echo " [setup] Checking T2T genome..."
if ! [ -f 2_source/T2T.fa ]; then
    echo "          - Downloading T2T genome [https://hgdownload.gi.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz]"
    wget -q -O 2_source/T2T.fa.gz https://hgdownload.gi.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
    gunzip -q 2_source/T2T.fa.gz
fi

# Download the 1000G T2T bcf
echo " [setup] Checking 1KG T2T bcf..."
if ! [ -f 1_plink/1KG.T2T.bcf.gz ]; then
    echo "          - Downloading 1KG T2T bcf...   [https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz]"
    wget -O 1_plink/1KG.T2T.bcf.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz
    wget -q -O 1_plink/1KG.T2T.bcf.gz.csi https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1/1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.3202.bcf.gz.csi
fi

# Download plink2
echo " [setup] Checking plink2..."
if ! [ -f 1_plink/plink2 ]; then
    echo "          - Downloading plink2...   [https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20250129.zip]"
    wget -q -O 1_plink/plink2.zip https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_avx2_20250129.zip
    unzip 1_plink/plink2.zip -d 1_plink
    rm 1_plink/plink2.zip 1_plink/vcf_subset 
fi

# Run plink to make pgen files
echo " [setup] Checking pgen files..."
if ! [ -f 1_plink/1KG.T2T.pvar ]; then
    echo "          - Running plink on 1KG to make pgen files... [see 1_plink/1KG.T2T.log]"
    echo
    1_plink/plink2 \
        --bcf 1_plink/1KG.T2T.bcf.gz \
        --geno 0.05 \
        --hwe  1e-10000 \
        --maf  0 \
        --make-pgen \
        --max-alleles 1000 \
        --new-id-max-allele-len 10000 missing \
        --not-chr chrX,chrY,chrM \
        --out 1_plink/1KG.T2T \
        --output-chr chrM \
        --set-all-var-ids '@_#' 
        # --silent
        echo
fi

# Download GWAS summary stats
echo " [setup] Checking GWAS summary files..."
if ! [ -f 2_source/GWAS.tsv ]; then
    echo "          - Downloading GWAS summary files...   [https://api.kpndataregistry.org/api/d/7j5797]"
    wget -q -O 2_source/GWAS_raw.zip https://api.kpndataregistry.org/api/d/7j5797
    unzip -q 2_source/GWAS_raw.zip -d 2_source/
    gunzip -q 2_source/GP2_euro_ancestry_meta_analysis_2024/GP2_ALL_EUR_ALL_DATASET_HG38_12162024.txt.gz -c > 2_source/GWAS_raw.tsv
    rm 2_source/GP2_euro_ancestry_meta_analysis_2024 -r
    rm -rf 2_source/__MACOSX/
    conda run -n PROX --no-capture-output python -c "
import numpy as np
import pandas as pd
GWAS = pd.read_csv('2_source/GWAS_raw.tsv', sep='\t')
GWAS['logp'] = -np.log(GWAS['p_value'])
GWAS['sig'] = GWAS['p_value'] < 5*10**-8
GWAS.to_csv('2_source/GWAS.tsv', sep='\t', index=None)"
fi

# Download GWAS leads and process
echo " [setup] Checking GWAS leads..."
if ! [ -f 2_source/GP2_leads.tsv ]; then
    echo "          - Downloading GWAS leads...   [https://www.medrxiv.org/content/medrxiv/early/2025/03/17/2025.03.14.24319455/DC2/embed/media-2.xlsx?download=true]"
    wget -q -O 2_source/GP2_leads.xlsx https://www.medrxiv.org/content/medrxiv/early/2025/03/17/2025.03.14.24319455/DC2/embed/media-2.xlsx?download=true
    conda run -n PROX --no-capture-output python -c "
import pysam
import pandas as pd
from pyliftover import LiftOver
from utils import liftover_variant
leads = pd.read_excel('2_source/GP2_leads.xlsx', sheet_name='Table S3')
fasta = pysam.FastaFile('2_source/T2T.fa')
lo = LiftOver('hg38','Hs1')
leads[['T2T','T2T_ref','T2T_alt']] = leads.apply(lambda r: liftover_variant('chr' + str(r['CHR']), r['BP'], r['Reference allele'], r['Alternative allele'], fasta, lo), result_type='expand', axis=1)
leads['ID'] = leads.apply(lambda r: 'chr' + str(r['CHR']) + '_' + str(r['T2T']), axis=1)
leads[['CHR','BP','Nearest Gene', 'Reference allele','Alternative allele','ID','Effect allele frequency']].rename(
        columns={'Reference allele':'ref','Alternative allele':'alt','Effect allele frequency':'AF'}).to_csv('2_source/GP2_leads.tsv', sep='\t', index=None)"
fi

# R-squared of leads to 1KG using plink
echo " [setup] Checking plink r-square output"
if ! [ -f 1_plink/linkage.vcor ]; then
    echo "          - Calculating R2 of leads to 1KG variants using plink..."
    # awk -F'\t' 'NR>1 { print $1 }' 2_source/GP2_leads.tsv | sed 's/:/_/g'  > 2_source/GP2_leads.id
    awk -F'\t' 'NR>1 { print $30 }' 2_source/GP2_leads.tsv > 2_source/GP2_leads.id
    echo
    1_plink/plink2 \
        --pfile 1_plink/1KG.T2T \
        --r2-phased \
        --ld-window-kb 2000 \
        --ld-window-r2 0.05 \
        --maf 0.001 \
        --max-alleles 2 \
        --ld-snp-list 2_source/GP2_leads.id \
        --out 1_plink/linkage
        # --silent
        echo
fi

# Check how many lead variants are in 1KG T2T after lifting to T2T
echo " [setup] Finding 1KG & GWAS lead overlap..."
if ! [ -f 3_variants/leads_in_1kg.txt ]; then
    grep -v '#' 1_plink/1KG.T2T.pvar | awk '{ print $3 }' | sort > 3_variants/all_hs1_ids.txt
    sort 2_source/GP2_leads.id > 3_variants/all_lead_ids.txt
    join 3_variants/all_lead_ids.txt 3_variants/all_hs1_ids.txt > 3_variants/leads_in_1kg.txt
    N=$( wc -l 3_variants/leads_in_1kg.txt)
    echo " [setup] Found $N leads in 1000 Genomes after lifting to T2T and converting with plink"
fi

# Check wether donor variants are copied here
echo " [setup] Checking donor bcf prescense..."
if ! [ -f 2_source/donor_vars.bcf ]; then
    echo "          - Copying donor variant bcf localy..."
    cp $DONOR_BCF 2_source/donor_vars.bcf
fi

echo
echo " Setup done."
echo
echo " ======================================================================================================================"
echo
