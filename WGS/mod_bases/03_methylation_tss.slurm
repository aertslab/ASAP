#!/bin/bash
#SBATCH -n 1
#SBATCH -c 15
#SBATCH --time=00:30:00
#SBATCH -A XX
#SBATCH -p XX
#SBATCH --cluster wice
#SBATCH --job-name=modkit_tss
#SBATCH --array=1-125%50
#SBATCH -o logs/modkit_%x_%j.out
#SBATCH -e logs/modkit_%x_%j.err

# script to calculate average TSS methylation for top and bottom 10% expressed genes 
# in a given cell type, based on a bedmethyl file and a TSS bed file.

# change the path to the singularity container/environment as needed. 
SING=/singularity_images/images_Seppe/pycistopic_2.0a.sif 

SAMPLES=(
ASA_001B
ASA_002B
ASA_004B
ASA_005B
ASA_011B
ASA_015B
ASA_016B
ASA_018B
ASA_020B
ASA_021B
ASA_022B
ASA_024B
ASA_025B
ASA_029B
ASA_031B
ASA_032B
ASA_033B
ASA_034B
ASA_035B
ASA_040B
ASA_041B
ASA_042B
ASA_043B
ASA_046B
ASA_048B
ASA_049B
ASA_050B
ASA_056B
ASA_057B
ASA_058B
ASA_059B
ASA_060B
ASA_064B
ASA_066B
ASA_070B
ASA_079B
ASA_080B
ASA_081B
ASA_082B
ASA_083B
ASA_084B
ASA_085B
ASA_088B
ASA_091B
ASA_092B
ASA_093B
ASA_094B
ASA_095B
ASA_105B
ASA_107B
ASA_108B
ASA_109B
ASA_111B
ASA_112B
ASA_113B
ASA_114B
ASA_115B
ASA_116B
ASA_117B
ASA_118B
ASA_119B
ASA_120B
ASA_121B
ASA_122B
ASA_123B
ASA_124B
ASA_125B
ASA_126B
ASA_127B
ASA_128B
ASA_129B
ASA_130B
ASA_131B
ASA_132B
ASA_133B
ASA_134B
ASA_135B
ASA_136B
ASA_138B
ASA_139B
ASA_140B
ASA_141B
ASA_142B
ASA_143B
ASA_144B
ASA_145B
ASA_148B
ASA_149B
ASA_150B
ASA_153B
ASA_154B
ASA_155B
ASA_158B
ASA_159B
ASA_160B
ASA_161B
ASA_162B
ASA_163B
ASA_164B
ASA_165B
ASA_166B
ASA_167B
ASA_168B
ASA_169B
ASA_170B
ASA_171B
ASA_172B
ASA_173B
ASA_174B
ASA_175B
ASA_176B
ASA_177B
ASA_178B
ASA_179B
ASA_180B
ASA_181B
ASA_182B
ASA_183B
ASA_184B
ASA_185B
ASA_186B
ASA_187B
ASA_188B
ASA_189B
ASA_190B
)

sample=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

# try input file - it should be the bed files created by 01_modkit_bed.slurm
if [ -f /mod_bases/${sample}_cpg.bed.gz ]; then
    echo "Input file found: /mod_bases/${sample}_cpg.bed.gz"
    file="/mod_bases/${sample}_cpg.bed.gz"
else
    echo "Input file not found for sample ${sample}. Skipping."
    exit 1
fi

# potentially add full path to the 03_methylation_tss_mean.py script if needed
apptainer run -B /staging,/lustre1,/data,/local_scratch,/scratch $SING python 03_methylation_tss_mean.py --sample $file --cell_type Oligo

