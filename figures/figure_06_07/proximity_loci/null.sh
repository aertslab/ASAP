#!/usr/bin/env bash
set -euo pipefail

R2_THRESHOLD=""
PASSTHROUGH_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --r2)
            R2_THRESHOLD="$2"
            shift 2
            ;;
        --peaks-df)
            PEAKS_DF="$2"
            shift 2
            ;;
        *)
            PASSTHROUGH_ARGS+=("$1")
            shift
            ;;
    esac
done

echo
echo " Running variant score enrichment analysis"
echo
echo " I."
echo " ---------------------------------------------------------------------------------------------------------"
echo " Derive MAF and CHR null distribution of lead variants"
echo " Oversample variants from the genome that create the same distribution"
echo

# Prep variant pool for bootstrapping null enrichment leads
echo " [null] Prepping variant pool for bootstrapping null enrichment leads... [1_plink/1KG.T2T.tsv]"
if ! [ -f 1_plink/1KG.T2T.tsv ]; then
    echo -e 'ID\tCHR\tBP\tMAF' > 1_plink/1KG.T2T.tsv
    grep -v '#' 1_plink/1KG.T2T.pvar | awk -F'\t' '{ print $1"_"$2";"$1";"$2";"$6 }' | awk -F';' '{ print $1"\t"$2"\t"$3"\t"$6 }' | sed 's/MAF=//g' >> 1_plink/1KG.T2T.tsv
fi

#
echo " [null] Getting matched MAF and CHR distribution of lead variants... [6_enrichment/null_variants.tsv]"
if ! [ -f 6_enrichment/null_variants.tsv ]; then
    conda run -n PROX --no-capture-output python null_variants.py --sample 1000
fi

# R-squared of null leads to 1KG using plink
echo " [null] R-squared of null leads to 1KG using plink... [1_plink/null.vcor]"
if ! [ -f 1_plink/null.vcor ] || ! [ -f 6_enrichment/null_variants.id ]; then
    echo "          - Calculating R2 of leads to 1KG variants using plink..."
    echo
    awk -F'\t' 'NR>1 { print $4 }' 6_enrichment/null_variants.tsv > 6_enrichment/null_variants.id
    1_plink/plink2 \
        --pfile 1_plink/1KG.T2T \
        --r2-phased \
        --ld-window-kb 2000 \
        --ld-window-r2 0.05 \
        --maf 0.001 \
        --max-alleles 2 \
        --ld-snp-list 6_enrichment/null_variants.id \
        --out 1_plink/null
    echo
fi

# Check how many null variants are in 1KG T2T after lifting to T2T
echo " [null] Finding 1KG & GWAS null lead overlap... [6_enrichment/nulls_in_1kg.txt]"
if ! [ -f 6_enrichment/nulls_in_1kg.txt ]; then
    sort 6_enrichment/null_variants.id > 6_enrichment/all_null_ids.txt
    join 6_enrichment/all_null_ids.txt 3_variants/all_hs1_ids.txt > 6_enrichment/nulls_in_1kg.txt
    N=$( wc -l 6_enrichment/nulls_in_1kg.txt)
    echo " [null] Found $N null leads in 1000 Genomes"
fi


echo
echo " II."
echo " ---------------------------------------------------------------------------------------------------------"
echo " Create loci around sampled variants"
echo " using the same algorithm as used on the leads"
echo

# Run proximity block generation from earlier but for null variants
echo " [null] Proximity block generation for null variants... [6_enrichment/regions.tsv]"
if ! [ -f 6_enrichment/regions.tsv ]; then

    CMD=(conda run -n PROX --no-capture-output python prox.py "${PASSTHROUGH_ARGS[@]}")
    if [[ -n "$R2_THRESHOLD" ]]; then
        CMD+=(--r2-threshold "$R2_THRESHOLD")
    CMD+=(--leads "6_enrichment/null_variants.tsv")
    CMD+=(--r2 "1_plink/null.vcor")
    CMD+=(--out-dir "6_enrichment")
    CMD+=(--variants-dir "6_enrichment")
    CMD+=(--variants-in-1kg "6_enrichment/nulls_in_1kg.txt")
    fi

    "${CMD[@]}"
fi


echo
echo " III."
echo " ---------------------------------------------------------------------------------------------------------"
echo " Get statistics for each lead locus and null locus:"
echo " SNP density, size, mean MAF"
echo " Save bootstrap null loci sets matching lead loci statistics distributions"
echo


# Creating bed of lead loci
echo " [stats] Creating bed of lead loci [6_enrichment/lead_regions.bed]"
awk 'NR>1 { print $1"\t"$11"\t"$12"\t"$2 }' 4_regions/regions.tsv | bedtools sort -i > 6_enrichment/lead_regions.bed

# Creating bed of null loci
echo " [stats] Creating bed of null loci [6_enrichment/null_regions.bed]"
awk 'NR>1 { print $1"\t"$11"\t"$12"\t"$2 }' 6_enrichment/regions.tsv | bedtools sort -i > 6_enrichment/null_regions.bed

echo " [stats] Creating 1kg variant bed file from bcf... [1_plink/1KG.T2T.bed]"
if ! [ -f 1_plink/1KG.T2T.bed ]; then
    grep -v '#' 1_plink/1KG.T2T.pvar | awk '{ print $1";"$2";"$2+1";"$3";"$6 }' | awk -F';' '{ print $1"\t"$2"\t"$3"\t"$4"\t"$7 }'  | sed 's/MAF=//g' | bedtools sort -i > 1_plink/1KG.T2T.bed
    #bcftools view 1_plink/1KG.T2T.bcf.gz | grep -v '#' | awk -F'\t' '{ print $1";"$2";"$2+1";"$1"_"$2";"$8 }' | awk -F';' '{ print $1"\t"$2"\t"$3"\t"$4"\t"$7 }'  | sed 's/MAF=//g' > 1_plink/1KG.T2T.bed
fi

echo " [stats] Intersecting lead loci with variant bed using bedtools... [6_enrichment/lead_loci_SNP.tsv]"
if ! [ -f 6_enrichment/lead_loci_SNP.tsv ]; then
    echo -e "chr\tvariant\tvariant_id\tstart\tstop\tlocus_id\tMAF" > 6_enrichment/lead_loci_SNP.tsv
    bedtools intersect -b 6_enrichment/lead_regions.bed -a 1_plink/1KG.T2T.bed -wa -wb | awk '{ print $1"\t"$2"\t"$4"\t"$7"\t"$8"\t"$9"\t"$5 }' >> 6_enrichment/lead_loci_SNP.tsv
fi

echo " [stats] Intersecting null loci with variant bed using bedtools... [6_enrichment/null_loci_SNP.tsv]"
if ! [ -f 6_enrichment/null_loci_SNP.tsv ]; then
    echo -e "chr\tvariant\tvariant_id\tstart\tstop\tlocus_id\tMAF" > 6_enrichment/null_loci_SNP.tsv
    bedtools intersect -b 6_enrichment/null_regions.bed -a 1_plink/1KG.T2T.bed -wa -wb | awk '{ print $1"\t"$2"\t"$4"\t"$7"\t"$8"\t"$9"\t"$5 }' >> 6_enrichment/null_loci_SNP.tsv
fi

echo
echo " [enrich] Get null-loci lead variants... [6_enrichment/draw_dict.pkl]"
if ! [ -f 6_enrichment/draw_dict.pkl ]; then
    conda run -n PROX --no-capture-output python null_loci.py --peaks-df "$PEAKS_DF"
fi


echo
echo " Loci bootstrap sets construction done."
echo
echo " ======================================================================================================================"
echo
