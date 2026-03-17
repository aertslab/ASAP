# WGS preprocessing

- for the general processing of raw sequencing data to aligned CRAM file, please refer to the following GitHub repository: [NanoWGS](https://github.com/AlexanRNA/nanowgs)
- in general, if R or Python environment is used, the environment info is summarised with `pip list` command in each file

## Modified bases extraction
- to extract modified bases from CRAM files in `./mod_bases`
- modkit version was dcontainerised and shared via dockerhub: `alexanrna/modkit:v0.4.3`
- contains scripts to overlap methylation info with TSS and CTCF sites

## meQTL analysis
- to rerun the analysis, follow the code `./meQTL` directory, the steps are ordered numerically. GPU is recommended for QTL mapping (`06_qtl_mapping.slurm`)
- tensorQTL v1.0.10 was used for mapping