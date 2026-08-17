# Load bcftools and bedtools (depends on system)
module load BCFtools/1.22-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0

# resource files ands plink
./setup.sh --donor-bcf /staging/leuven/stg_00090/ASA/analysis/20241004_patient_demultiplexing/WGS_chm13_BCFs_merged/WGS_chm13_BCFs.missing_to_ref.norm.remove_overlaps.bcf

# Create loci around lead variants
./prox.sh --r2 0.6

echo
echo " Enrichment analysis:"
echo

# Create null-matched loci
./null.sh --r2 0.6 --peaks-df /lustre1/project/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/2026_full_dataset/all/out/combined_consensus_peaks_500bp.bed

# Run actual enrichtment analysis
conda run -n PROX --no-capture-output python null_tail.py \
                  --score-dataframe /staging/leuven/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/2026_full_dataset/20260526_CREsted_full/20260719_CREsted_variant_scoring/cutsite_model_finetuned11/variant_scores_long.parquet \
                  --metric-column logfc \
                  --threshold 0.5

conda run -n PROX --no-capture-output python null_peaks.py
