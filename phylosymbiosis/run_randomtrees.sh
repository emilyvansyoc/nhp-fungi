#!/bin/bash

##### run Andrew's script to generate random trees
echo "###### GENERATING RANDOM TREES ####"
# activate conda environment
eval "$(conda shell.bash hook)"
conda activate py2

## ---- set variables ----

# set working directory
DIR=phylosymbiosis/

# number of random trees to generate
NRAND=100000

# tree labels (comma separated list with no spaces)
LABS=$1

# name for random tree (outside variable to be easy to run)
RAND=$2

## ---- run random trees ----

# make temporary directory for random trees
TMP=$DIR/TEMP
mkdir -p $TMP

# run random trees
python $DIR/tree_random.py -n $NRAND -i $LABS -o $TMP

# concatenate trees into one file
for i in $TMP/*
do echo $i
    cat $i >> $RAND
done

# remove temp directory
rm -rf $TMP

# print finished message
echo "done making random trees"