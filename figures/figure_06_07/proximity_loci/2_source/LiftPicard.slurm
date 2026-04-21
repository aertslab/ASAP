#!/bin/bash -l
#SBATCH --clusters="wice"
#SBATCH --account="lp_big_wice_cpu"
#SBATCH --job-name="Lift GWAS Picard"
#SBATCH --partition="dedicated_big_bigmem"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --time=2:00:00
#SBATCH --gpus-per-node=0
#SBATCH -o "/staging/leuven/stg_00090/ASA/analysis/analysis_bram/Sigalova_Pancikova_DeMan_2026/95.1.fig67/2_source/picard.log"

cd /staging/leuven/stg_00090/ASA/analysis/analysis_bram/Sigalova_Pancikova_DeMan_2026/95.1.fig67/2_source

module load BCFtools/1.22-GCC-10.3.0

### Correct chr23/24 -> X,Y, convert .tsv to .vcf, fix vcf header, noramlise, index
awk 'NR==1{print; next}
{
  if ($1==23) $1="chrX";
  else if ($1==24) $1="chrY";
  else $1="chr"$1;
  print
}' OFS="\t" GWAS.tsv | awk 'BEGIN {
    OFS="\t";
    print "##fileformat=VCFv4.2";
    print "##reference=hg38";
    print "##INFO=<ID=BETA,Number=1,Type=Float,Description=\"Effect size\">";
    print "##INFO=<ID=PVAL,Number=1,Type=Float,Description=\"P-value\">";
    print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO";
}
NR>1 {
    chrom=$1;
    pos=$2;
    id=$3;
    alt=$4;
    ref=$5;
    beta=$9;
    pval=$8;

    info="BETA="beta";PVAL="pval;

    print chrom, pos, id, ref, alt, ".", "PASS", info;
}' | bcftools reheader -f /staging/leuven/res_00001/genomes/homo_sapiens/hg38_ucsc/fasta/hg38.fa.fai | bcftools norm \
    -f /staging/leuven/res_00001/genomes/homo_sapiens/hg38_ucsc/fasta/hg38.fa \
    -m -any \
    -Oz -o GWAS.hg38.gz

bcftools index GWAS.hg38.gz

module load picard/3.1.1-Java-19.0.2

liftover_vcf_from_hg38_to_chm13 () {
    local vcf_input_filename="${1}";
    local vcf_output_filename="${2}";

    local hg38_to_chm13_chroms_chain_filename='/staging/leuven/res_00001/genomes/homo_sapiens/CHM13v2_maskedY_rCRS/chain_files/hg38ToHs1.over.chain';
    local chm13_reference_fasta_filename='/staging/leuven/res_00001/genomes/homo_sapiens/CHM13v2_maskedY_rCRS/fasta/chm13v2.0_maskedY_rCRS.fa';

    java \
	-Djava.io.tmpdir="${TMPDIR:-/tmp}" \
        -jar ${EBROOTPICARD}/picard.jar \
        LiftoverVcf \
            I="${vcf_input_filename}" \
            O="${vcf_output_filename}" \
            CHAIN="${hg38_to_chm13_chroms_chain_filename}" \
            REJECT="${vcf_output_filename%.vcf}.rejected.vcf.gz" \
            WARN_ON_MISSING_CONTIG=true \
            LOG_FAILED_INTERVALS=true \
            WRITE_ORIGINAL_POSITION=true \
            WRITE_ORIGINAL_ALLELES=true \
            LIFTOVER_MIN_MATCH=0.9 \
            ALLOW_MISSING_FIELDS_IN_HEADER=false \
            RECOVER_SWAPPED_REF_ALT=true \
            REFERENCE_SEQUENCE="${chm13_reference_fasta_filename}"

    return $?;
}

liftover_vcf_from_hg38_to_chm13 GWAS.hg38.gz GWAS.hs1.gz

gunzip GWAS.hs1.gz

grep -v '#' GWAS.hs1 | awk '{ print $1"_"$2"_"$4"_"$5"\t"$1"\t"$2"\t"$4"\t"$5 }' | awk -F'\t' '{ print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5 }' > GWAS.hs1.tsv
