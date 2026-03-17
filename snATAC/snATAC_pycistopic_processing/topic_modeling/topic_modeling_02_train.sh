#!/bin/bash -l
brain_region=$1
nt=$2

temp=/scratch/$brain_region/topic_modelling_$nt/
corpus_folder=/scratch/$brain_region/344/corpus/
mkdir $temp

export MALLET_MEMORY="400g"


/data/mallet/Mallet-202108/bin/mallet train-topics \
    --input $corpus_folder/corpus \
    --num-topics $nt \
    --alpha 50 \
    --beta 0.1 \
    --optimize-interval 0 \
    --num-threads 30 \
    --output-doc-topics $temp/$nt.topics_doctopics.txt \
    --word-topic-counts-file $temp/$nt.word_topics_counts.txt \
    --output-topic-keys $temp/$nt.topics__topickeys.txt \
    --num-iterations 200 \
    --doc-topics-threshold 0.0 \
    --random-seed 555