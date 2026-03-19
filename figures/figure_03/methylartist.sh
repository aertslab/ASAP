#/bin/bash

mamba activate methylartist
methylartist locus -b ASA_021B_phased.cram -m m --motif CG --ref /reference/chm13_v2.0_maskerY.rCRS/fasta/chm13v2.0_maskedY_rCRS.fa \
-i chr1:209125800-209127630 --primary_only --outfile SYT14_ASA_021B.var.pdf --variants chr1_SYT14_var.vcf.gz --splitvar SYT1 \
--samplepalette magma -l chr1:209126280-209126780  --variantsize 5 --skip_raw_plot --highlightpalette Pastel1