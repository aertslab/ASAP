#!/bin/bash -l
covs=$2
var=$1
peaks=$3
br=$4

SIF=/software/tensorqtl_1.0.10.sif

singularity exec --nv \
  -B /staging,/lustre1,/data,/vsc-hard-mounts,/local_scratch,/scratch,/user,/apps \
  $SIF python3 --version

if [ "$br" == "SN" ]; then
    cell_types=("Astro" "DopaN" "Endo" "GabaN" "Micro-PVM" "OPC" "Oligo" "neuron")
elif [ "$br" == "CC" ]; then
    cell_types=("Astro" "L23_IT" "Endo" "L4_IT" "L56_NP" "L5_ET" "L5_IT_A" "L5_IT_B" "L6_CT" "L6_IT_Car3" "L6_IT" "L6b" "Lamp5" "Micro-PVM" "OPC" "Oligo" "Pvalb" "Sncg" "Sst" "Vip")
else
    echo "Unknown brain region: $br"
fi

export var covs peaks br

# Run in parallel
parallel --jobs 3 "echo Running job for cell type: {} && \
    singularity exec --nv \
      -B /staging,/lustre1,/data,/vsc-hard-mounts,/local_scratch,/scratch,/user,/apps \
      $SIF python3 tensorQTL_parquet.py \
      -br \"$br\" -ct {} -vt \"$var\" -cov \"$covs\" -peaks \"$peaks\"" ::: "${cell_types[@]}"