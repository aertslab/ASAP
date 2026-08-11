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
        *)
            PASSTHROUGH_ARGS+=("$1")
            shift
            ;;
    esac
done

CMD=(conda run -n PROX --no-capture-output python prox.py "${PASSTHROUGH_ARGS[@]}")
if [[ -n "$R2_THRESHOLD" ]]; then
    CMD+=(--r2-threshold "$R2_THRESHOLD")
fi

echo

"${CMD[@]}"



# Convert donor variant bcf to tsv
echo " [overlap] Converting donor variant bcf to tsv... [creating 3_variants/donor_vars.tsv]"
if ! [ -f 3_variants/donor_vars.tsv ]; then
    module load BCFtools/1.22-GCC-10.3.0
    bcftools view 2_source/donor_vars.bcf | grep -v '#' | awk '{ print $1"\t"$2"\t"$4"\t"$5"\t"$8 }' | awk -F'\t' '{ print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5 }' > 3_variants/donor_vars.tsv
fi

# Find variants in loci using bedtools
echo " [overlap] Find variants in loci using bedtools... [creating 3_variants/donor_subset.tsv]"
if ! [ -f 3_variants/donor_subset.tsv ]; then
    module load BEDTools/2.30.0-GCC-10.3.0
    echo -e "chr\tbp\tref\talt\tlead_id" > 3_variants/donor_subset.tsv
    bedtools intersect -b <(awk 'NR>1 {print $0}' 4_regions/regions.tsv) -a 3_variants/donor_vars.tsv -wa -wb | awk '{ print $1"\t"$2"\t"$4"\t"$5"\t"$10 }' >> 3_variants/donor_subset.tsv
fi

echo
echo "Done."
echo