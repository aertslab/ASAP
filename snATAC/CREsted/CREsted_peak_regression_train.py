import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ad
import tensorflow as tf
import wandb
import keras
import crested
import pickle 
import sys

analysis_dir = "/staging/leuven/stg_00090/ASA/analysis/"

# read adata path and model name from command line

adata_path = sys.argv[1]
model_name = sys.argv[2]

#adata_path = f"{analysis_dir}/analysis_Olga/3_T2T_analysis/2_CREsted_models/peak_regression/adata_cc_norm.h5ad"
#model_name = "deepPeak_CC_basemodel"

############################


# Set the genome
genome = crested.Genome(
    "/staging/leuven/res_00001/genomes/homo_sapiens/CHM13v2_maskedY_rCRS/fasta/chm13v2.0_maskedY_rCRS.fa", 
    "/staging/leuven/res_00001/genomes/homo_sapiens/CHM13v2_maskedY_rCRS/fasta/chm13v2.0_maskedY_rCRS.chrom.sizes"
)
crested.register_genome(
    genome
)

# Load adata
adata = ad.read_h5ad(adata_path)


# Set data module
datamodule = crested.tl.data.AnnDataModule(
    adata,
    batch_size=256,  # lower this if you encounter OOM errors
    max_stochastic_shift=3,  # optional data augmentation to slightly reduce overfitting
    always_reverse_complement=True,  # default True. Will double the effective size of the training dataset.
)


# Load chrombpnet-like architecture for a dataset with 2114bp regions and CC cell types
model_architecture = crested.tl.zoo.dilated_cnn(
    seq_len=2114, num_classes=len(list(adata.obs_names))
)


# Create configuration for peak regression with a weighted cosine mse log loss function
optimizer = keras.optimizers.Adam(learning_rate=1e-3)
loss = crested.tl.losses.CosineMSELogLoss(max_weight=100) # for predicting counts set multiplier=1

metrics = [
    keras.metrics.MeanAbsoluteError(),
    keras.metrics.MeanSquaredError(),
    keras.metrics.CosineSimilarity(axis=1),
    crested.tl.metrics.PearsonCorrelation(),
    crested.tl.metrics.ConcordanceCorrelationCoefficient(),
    crested.tl.metrics.PearsonCorrelationLog(),
    crested.tl.metrics.ZeroPenaltyMetric(),
]

alternative_config = crested.tl.TaskConfig(optimizer, loss, metrics)


# Setup the trainer
trainer = crested.tl.Crested(
    data=datamodule,
    model=model_architecture,
    config=alternative_config,
    project_name="deepPeak",  # change to your liking
    run_name=model_name,  # change to your liking
    logger="wandb",  # or None, 'dvc', 'tensorboard'
    seed=7,  # For reproducibility
)

# Train the model
trainer.fit(
    epochs=60,
    learning_rate_reduce_patience=3,
    early_stopping_patience=6,
)