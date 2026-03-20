# snATAC-seq analysis
- preprocessing is under `atac_preprocessing/` (includes barcode trimming, mapping and merging bams per sample)
- pycisTopic analysis per brain region is under `snATAC_pycistopic_processing/` (includes: QC, topic modeling, annotation, pseudobulking and peak calling, DARs, creation of matrices for caQTLs) 
- caQTL analysis per brain region `caQTLs/` (includes: matrix preparation, covariat preparation, tensorQTL)
- ASCA-related analysis is under `wasp_pipeline/` (under construction)