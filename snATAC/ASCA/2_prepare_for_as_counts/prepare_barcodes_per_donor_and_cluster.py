import os
import glob
import pysam
import pandas as pd

# Load and label CC data
path_cc = "/lustre1/project/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/CC/out/CC_cell2cluster_20250121.tsv"
df1 = pd.read_csv(path_cc)
df1["brain_region"] = "CC"

# Load and label SN data
path_sn = "/lustre1/project/stg_00090/ASA/analysis/2024_T2T_ATAC_analysis/SN/out/SN_cell2cluster_20241009.tsv"
df2 = pd.read_csv(path_sn)
df2["brain_region"] = "SN"

# Combine data
df = pd.concat([df1, df2])
df.rename({"cell_barcode_input": "cell_barcode"}, axis=1, inplace=True)

# Extract donor ID from cluster field
df["donor_id"] = df["cluster"].str.extract(r'(ASA_[0-9]+)', expand=False)

# Extract cell type and clean up
df["cell_type"] = df["cluster"].str.split("_").str[0]
df = df[(df.cell_type != "unknown") & (df.cell_type != "doublet")]
df["cell_type"] = df.cell_type.str.replace("/", "").str.replace(" ", "_")

# Extract number of cells
df["n_cells"] = df["cluster"].str.extract(r'_([0-9]+)$', expand=False)

# Compose cluster name
df["cluster_name"] = df.brain_region + "." + df.cell_type + "." + df.donor_id + "." + df.n_cells

# Select final columns
df_sel = df[["cell_barcode", "cluster_name", "brain_region", "donor_id", "cell_type", "n_cells"]]

# Save output
outdir = "/lustre1/project/stg_00090/ASA/analysis/analysis_Olga/3_T2T_analysis/data/barcodes_per_donor"
os.makedirs(outdir, exist_ok=True)
df_sel.to_csv(os.path.join(outdir, "selected_barcodes_combined.tsv"), sep="\t", index=False)
