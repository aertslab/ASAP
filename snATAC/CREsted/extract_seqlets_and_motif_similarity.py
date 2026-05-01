import tfmindi as tm
import os
import re
import numpy as np
from tqdm import tqdm
import polars as pl
import pickle

############################
### Load motifs
############################

motif_collection_dir = tm.fetch_motif_collection()
motif_annotations_file = tm.fetch_motif_annotations()

# load sampled motifs
motif_samples_path = "sampled_motifs.txt"
with open(motif_samples_path) as f:
    motif_names = [line.strip() for line in f.readlines()]

print(f"Loaded {len(motif_names)} motifs")

motif_collection = tm.load_motif_collection(
    motif_collection_dir,
    motif_names=motif_names
)
motif_annotations = tm.load_motif_annotations(motif_annotations_file)
motif_to_db = tm.load_motif_to_dbd(motif_annotations)

outdir = "seqlets_v1_2"
os.makedirs(outdir, exist_ok=True)

############################
### Process seqlets
############################

brain_regions = ["CC", "SN"]
alleles = ["high_allele", "low_allele"]

for brain_region in brain_regions:

    print(f"\n=== Processing {brain_region} (high + low combined) ===")

    VAR_FOLDER = f"contribution_scores/{brain_region}/variant_lists"
    CONTRIB_FOLDERS = {
        "high_allele": f"contribution_scores/{brain_region}/high_allele",
        "low_allele": f"contribution_scores/{brain_region}/low_allele",
    }

    out_h5ad = f"{outdir}/seqlets_{brain_region}_combined.h5ad"
    if os.path.exists(out_h5ad):
        print(f"[SKIP] {out_h5ad} already exists")
        continue

    ############################
    # Load variant annotations
    ############################

    print("Loading variant annotations...")

    var_list = []

    class_names = [
        re.match(r"(.+?)_oh\.npz$", f).group(1)
        for f in os.listdir(CONTRIB_FOLDERS["high_allele"])
        if f.endswith("_oh.npz")
    ]

    for c in tqdm(class_names):
        var_list.append(
            pl.read_csv(os.path.join(VAR_FOLDER, f"{c}.csv"))
            .with_columns(pl.lit(c).alias("cell_type"))
        )

    var_df = pl.concat(var_list)

    ############################
    # Load contrib + OH (both alleles)
    ############################

    contrib_all = []
    oh_all = []

    print("Loading contribution scores (both alleles)...")

    for allele in alleles:
        folder = CONTRIB_FOLDERS[allele]

        for c in tqdm(class_names, desc=allele):
            contrib = np.load(os.path.join(folder, f"{c}_contrib.npz"))["arr_0"]
            oh = np.load(os.path.join(folder, f"{c}_oh.npz"))["arr_0"]

            contrib_all.append(contrib)
            oh_all.append(oh)

    oh = np.concatenate(oh_all)
    contrib = np.concatenate(contrib_all)

    # transpose to (n_examples, width, 4)
    oh = oh.transpose(0, 2, 1).astype("float32")
    contrib = contrib.transpose(0, 2, 1).astype("float32")

    print("  oh shape:", oh.shape)
    print("  contrib shape:", contrib.shape)

    ############################
    # Expand variant table to match allele stacking
    ############################

    var_df = pl.concat(
        [
            var_df.with_columns(pl.lit(a).alias("allele"))
            for a in alleles
        ],
        how="vertical"
    ).with_row_index("example_idx")

    # --- Sanity check ---
    assert var_df.height == oh.shape[0], (
        f"Mismatch: var_df={var_df.height}, oh={oh.shape[0]}"
    )

    ############################
    # Extract seqlets
    ############################

    print("Extracting seqlets (combined)...")

    seqlets_df, seqlets_matrices = tm.pp.extract_seqlets(
        contrib=contrib,
        oh=oh,
        threshold=0.01,
        additional_flanks=3
    )

    seqlets_df = pl.from_pandas(seqlets_df)

    print(f"Extracted {len(seqlets_matrices)} seqlets")

    ############################
    # Annotate seqlets (CORRECT WAY)
    ############################

    seqlets_df = (
        seqlets_df
        .join(var_df, on="example_idx", how="left")
        .with_columns(
            pl.len()
            .over("peak_id", "cell_type", "allele")
            .alias("n_seqlets_per_peak")
        )
    )

    # --- Sanity checks ---
    assert seqlets_df["cell_type"].null_count() == 0
    assert seqlets_df["allele"].null_count() == 0
    assert seqlets_df["example_idx"].null_count() == 0

    print("Seqlet annotation sanity checks passed")

    seqlets_df.write_csv(
        f"{outdir}/seqlets_{brain_region}_combined_annotated.csv"
    )

    seqlets_df_pd = seqlets_df.to_pandas()

    with open(
        f"{outdir}/seqlet_matrices_{brain_region}_combined.pkl", "wb"
    ) as f:
        pickle.dump(seqlets_matrices, f, protocol=pickle.HIGHEST_PROTOCOL)

    ############################
    # Motif similarity
    ############################

    print("Calculating motif similarities...")

    sim_matrix = tm.pp.calculate_motif_similarity(
        seqlets_matrices,
        motif_collection,
        chunk_size=50000
    )

    print(f"Similarity matrix shape: {sim_matrix.shape}")

    adata = tm.pp.create_seqlet_adata(
        sim_matrix,
        seqlets_df_pd,
        seqlet_matrices=seqlets_matrices,
        oh_sequences=oh,
        contrib_scores=contrib,
        motif_collection=motif_collection,
        motif_annotations=motif_annotations,
        motif_to_dbd=motif_to_db,
    )

    tm.save_h5ad(adata, out_h5ad)

    print(f"Saved {out_h5ad}")
