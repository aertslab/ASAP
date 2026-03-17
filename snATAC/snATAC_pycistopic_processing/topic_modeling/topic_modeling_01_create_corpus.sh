#!/bin/bash -l
brain_region=$1
cto=$2
mallet=/data/mallet/Mallet-202108/bin/mallet
tmp_dir=/scratch/$brain_region/
outfile=$tmp_dir/corpus

mkdir $tmp_dir

python topic_modeling_01_create_corpus.py -cto $cto -m $mallet -out $outfile