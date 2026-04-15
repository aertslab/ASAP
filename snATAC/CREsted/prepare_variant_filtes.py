import os
import polars as pl

analysis_dir = "/staging/leuven/stg_00090/ASA/analysis/"

path = f"{analysis_dir}/2024_T2T_ATAC_analysis//20250518_combined_variant_scoring/ISM_ASCA_caQTL_combined_signif_with_eQTL_slopes.parquet"
df = pl.read_parquet(path)

outdir = f"{analysis_dir}/2024_T2T_ATAC_analysis/20250610_CREsted_variant_scoring_benchmark"
os.makedirs(outdir, exist_ok=True)

for brain_region in ["CC", "SN"]:
    df_sel = (
        df.filter(pl.col("brain_region") == brain_region)
        .with_columns(
            pl.col("variant_id")
            .str.split("_")
            .list[1]
            .cast(pl.Int64)
            .alias("variant_pos")
        )
        .select(
            "chr",
            "REGION.START",
            "REGION.END",
            "variant_pos",
            "TEST.SNP.REF.ALLELE",
            "TEST.SNP.ALT.ALLELE",
            "variant_id",
        )
        .unique()
    )

    print(f"Processing {brain_region} with {df_sel.shape[0]} variants")

    df_sel.write_csv(
        f"{outdir}/caQTL_ASCA_variants_{brain_region}.txt",
        separator="\t",
        include_header=False,
    )