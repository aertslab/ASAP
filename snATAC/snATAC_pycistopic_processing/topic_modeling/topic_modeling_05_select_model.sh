#!/bin/bash -l
br=$1
model_dir=/snATAC_analysis/$br/out/topic_modelling/models/

nt=$2
cto=$3
outdir=/snATAC_analysis/$br/out/topic_modelling/eval/$nt
mkdir $outdir

python topic_modeling_05_select_model.py \
    -i $cto \
    -tmd $model_dir \
    -m $nt \
    -o $outdir

mv $outdir/{$nt}_cto.pkl /snATAC_analysis/$br/out/{$br}_T2T_final_{$nt}topics_ct2.pkl