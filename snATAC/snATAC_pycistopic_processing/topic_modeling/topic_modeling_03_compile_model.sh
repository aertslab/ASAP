#!/bin/bash -l
brain_region=$1
nt=$2
cto=$3
tmp_dir=/scratch/$brain_region/topic_modelling_$nt/
out=/snATAC_analysis/$brain_region/out/topic_modelling/models/Topic{$nt}_test.pkl.gz

python topic_modeling_03_compile_model.py \
    -tmp $tmp_dir \
    -counts $tmp_dir/$nt.word_topics_counts.txt \
    -t $nt \
    -m /data/mallet/Mallet-202108/bin/mallet \
    -ncpu 10 \
    -save $out \
    -dt $tmp_dir/$nt.topics_doctopics.txt \
    -cto $cto \