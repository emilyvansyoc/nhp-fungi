#!/bin/bash

## run topological congruency script
# on QIIME2 dendrograms
# EVS 4/2025

# trap command+c to quit script and then 'set' to print commands
trap 'echo "Command+c, exiting script gracefully..."; exit 1' SIGINT
set -x

# location to script
SCRIPT="phylosymbiosis/run_topocong_script.R"

# set location to QIIME dendrograms that have been rooted in Figtree and rotated in R
INPUT_DIR="data/topocong_input/"

# set location to host tree (will be the same for all runs)
HOST_TREE="data/topocong_input/host_ITS2_rotated.newick"

# set main location for output
OUTPUT_DIR="data/topocong_output"

# set path to random trees
RAND_PATH="data/topocong_input/random_trees.tre"

# print help message and quit
Rscript $SCRIPT --help

########### OTU LEVEL ###########


# BRAY CURTIS
Rscript $SCRIPT -i $HOST_TREE -d $INPUT_DIR/bray_OTU_rotated_meanrare1837.newick -m false -r $RAND_PATH -o $OUTPUT_DIR/bray_OTU_meanrare1837_v2rotate/ 

# JACCARD
Rscript $SCRIPT -i $HOST_TREE -d $INPUT_DIR/jaccard_OTU_rotated_meanrare1837.newick -m false -r $RAND_PATH -o $OUTPUT_DIR/jaccard_OTU_meanrare1837 

########### GENUS LEVEL ###########

Rscript $SCRIPT -i $HOST_TREE -d $INPUT_DIR/bray_genusinR_rotated_meanrare1275.newick -m false -r $RAND_PATH -o $OUTPUT_DIR/bray_genusinR_meanrare1275 

Rscript $SCRIPT -i $HOST_TREE -d $INPUT_DIR/jaccard_genusinR_rotated_meanrare1275.newick -m false -r $RAND_PATH -o $OUTPUT_DIR/jaccard_genusinR_meanrare1275 
