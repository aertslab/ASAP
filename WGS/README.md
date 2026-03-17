# WGS preprocessing

- for the general processing of raw sequencing data to aligned CRAM file, please refer to the following GitHub repository: [NanoWGS](https://github.com/AlexanRNA/nanowgs)
- in general, if R or Python environment is used, the environment info is summarised with `pip list` command in each file

## mod_bases 
- to extract modified bases from CRAM files

## meQTL analysis
- to rerun the analysis, follow the code `./meQTL` directory, the steps are ordered numerically. GPU is recommended for QTL mapping (`06_qtl_mapping.slurm`)