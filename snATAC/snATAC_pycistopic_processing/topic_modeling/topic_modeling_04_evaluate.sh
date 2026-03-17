#!/bin/bash -l
brain_region=$1
model_dir=snATAC_analysis/$brain_region/out/topic_modelling/models/
nt=$2
cto=$3

outdir=/snATAC_analysis/$brain_region/out/topic_modelling/eval/$nt
mkdir $outdir

python topic_modeling_04_evaluate.py \
    -i $cto \
    -tmd $model_dir \
    -m $nt \
    -o $outdir