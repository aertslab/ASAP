#!/bin/bash
#SBATCH -n 1
#SBATCH -c 15
#SBATCH --mem 200G
#SBATCH --time=02:00:00
#SBATCH -A XX
#SBATCH -p XX
#SBATCH --cluster wice
#SBATCH --job-name=modbed_processing
#SBATCH --array=1-125%50
#SBATCH -o /logs/modbed_processing_%a_%j.out
#SBATCH -e /logs/modbed_processing_%a_%j.err
#SBATCH --mail-user=XX
#SBATCH --mail-type=FAIL,ARRAY_TASKS

# update paths to the input files as needed, procuded by ../modbases/01_modkit_bed.slurm 
SAMPLES=(./ASA_001B_cpg.bed.gz
./ASA_002B_cpg.bed.gz
./ASA_004B_cpg.bed.gz
./ASA_005B_cpg.bed.gz
./ASA_011B_cpg.bed.gz
./ASA_015B_cpg.bed.gz
./ASA_016B_cpg.bed.gz
./ASA_018B_cpg.bed.gz
./ASA_020B_cpg.bed.gz
./ASA_021B_cpg.bed.gz
./ASA_022B_cpg.bed.gz
./ASA_024B_cpg.bed.gz
./ASA_025B_cpg.bed.gz
./ASA_029B_cpg.bed.gz
./ASA_031B_cpg.bed.gz
./ASA_032B_cpg.bed.gz
./ASA_033B_cpg.bed.gz
./ASA_034B_cpg.bed.gz
./ASA_035B_cpg.bed.gz
./ASA_040B_cpg.bed.gz
./ASA_041B_cpg.bed.gz
./ASA_042B_cpg.bed.gz
./ASA_043B_cpg.bed.gz
./ASA_046B_cpg.bed.gz
./ASA_048B_cpg.bed.gz
./ASA_049B_cpg.bed.gz
./ASA_050B_cpg.bed.gz
./ASA_056B_cpg.bed.gz
./ASA_057B_cpg.bed.gz
./ASA_058B_cpg.bed.gz
./ASA_059B_cpg.bed.gz
./ASA_060B_cpg.bed.gz
./ASA_064B_cpg.bed.gz
./ASA_066B_cpg.bed.gz
./ASA_070B_cpg.bed.gz
./ASA_079B_cpg.bed.gz
./ASA_080B_cpg.bed.gz
./ASA_081B_cpg.bed.gz
./ASA_082B_cpg.bed.gz
./ASA_083B_cpg.bed.gz
./ASA_084B_cpg.bed.gz
./ASA_085B_cpg.bed.gz
./ASA_088B_cpg.bed.gz
./ASA_091B_cpg.bed.gz
./ASA_092B_cpg.bed.gz
./ASA_093B_cpg.bed.gz
./ASA_094B_cpg.bed.gz
./ASA_095B_cpg.bed.gz
./ASA_105B_cpg.bed.gz
./ASA_107B_cpg.bed.gz
./ASA_108B_cpg.bed.gz
./ASA_109B_cpg.bed.gz
./ASA_111B_cpg.bed.gz
./ASA_112B_cpg.bed.gz
./ASA_113B_cpg.bed.gz
./ASA_114B_cpg.bed.gz
./ASA_115B_cpg.bed.gz
./ASA_116B_cpg.bed.gz
./ASA_117B_cpg.bed.gz
./ASA_118B_cpg.bed.gz
./ASA_119B_cpg.bed.gz
./ASA_120B_cpg.bed.gz
./ASA_121B_cpg.bed.gz
./ASA_122B_cpg.bed.gz
./ASA_123B_cpg.bed.gz
./ASA_124B_cpg.bed.gz
./ASA_125B_cpg.bed.gz
./ASA_126B_cpg.bed.gz
./ASA_127B_cpg.bed.gz
./ASA_128B_cpg.bed.gz
./ASA_129B_cpg.bed.gz
./ASA_130B_cpg.bed.gz
./ASA_131B_cpg.bed.gz
./ASA_132B_cpg.bed.gz
./ASA_133B_cpg.bed.gz
./ASA_134B_cpg.bed.gz
./ASA_135B_cpg.bed.gz
./ASA_136B_cpg.bed.gz
./ASA_138B_cpg.bed.gz
./ASA_139B_cpg.bed.gz
./ASA_140B_cpg.bed.gz
./ASA_141B_cpg.bed.gz
./ASA_142B_cpg.bed.gz
./ASA_143B_cpg.bed.gz
./ASA_144B_cpg.bed.gz
./ASA_145B_cpg.bed.gz
./ASA_148B_cpg.bed.gz
./ASA_149B_cpg.bed.gz
./ASA_150B_cpg.bed.gz
./ASA_153B_cpg.bed.gz
./ASA_154B_cpg.bed.gz
./ASA_155B_cpg.bed.gz
./ASA_158B_cpg.bed.gz
./ASA_159B_cpg.bed.gz
./ASA_160B_cpg.bed.gz
./ASA_161B_cpg.bed.gz
./ASA_162B_cpg.bed.gz
./ASA_163B_cpg.bed.gz
./ASA_164B_cpg.bed.gz
./ASA_165B_cpg.bed.gz
./ASA_166B_cpg.bed.gz
./ASA_167B_cpg.bed.gz
./ASA_168B_cpg.bed.gz
./ASA_169B_cpg.bed.gz
./ASA_170B_cpg.bed.gz
./ASA_171B_cpg.bed.gz
./ASA_172B_cpg.bed.gz
./ASA_173B_cpg.bed.gz
./ASA_174B_cpg.bed.gz
./ASA_175B_cpg.bed.gz
./ASA_176B_cpg.bed.gz
./ASA_177B_cpg.bed.gz
./ASA_178B_cpg.bed.gz
./ASA_179B_cpg.bed.gz
./ASA_180B_cpg.bed.gz
./ASA_181B_cpg.bed.gz
./ASA_182B_cpg.bed.gz
./ASA_183B_cpg.bed.gz
./ASA_184B_cpg.bed.gz
./ASA_185B_cpg.bed.gz
./ASA_186B_cpg.bed.gz
./ASA_187B_cpg.bed.gz
./ASA_188B_cpg.bed.gz
./ASA_189B_cpg.bed.gz
./ASA_190B_cpg.bed.gz)

# update paths to singularity image and script if needed
SING=../software/singularity_images/images_Seppe/pycistopic_2.0a.sif
SCRIPT=02_modbed_processing.py

sample=${SAMPLES[$SLURM_ARRAY_TASK_ID-1]}

echo "Processing sample $sample"
apptainer run -B /staging/,/data/leuven/ \
$SING $SCRIPT --input_file $sample

