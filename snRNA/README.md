# snRNA-seq analysis
- all analysis is in `01_analysis`
- we first QC data from cellranger (`01_processing_10x`) and parsebiosciences (`01_processing_parsebio`)
- we then combine the data and finalyse annotations in `02_combination_10x_pb`
- to perform differential gene expression, go to `03_dge`
- to perform eQTL analysis, go to `04_QTL`
- there is also `plink` directory with plink scripts to prepare input for eQTL analysis

> **Note** <br>
The environmentds are either specified as `pip list` at the end of the script, or the docker container is mentioned at the beginning