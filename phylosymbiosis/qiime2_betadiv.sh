#!/bin/bash

# create directories for input and output files
INDIR="data/for-qiime2"
OUTPUT_DIR="data/qiime2-output/"
mkdir -p $OUTPUT_DIR

# First, run the R script to convert phyloseq to QIIME2 format
#Rscript phyloseq_to_qiime2.R

# Install QIIME2 using conda

# Download the QIIME2 environment file for the latest available version
cd $INDIR
#CONDA_SUBDIR=osx-64 conda env create -n qiime2-amplicon-2024.10 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-conda.yml



# Activate the environment
# Initialize conda for shell interaction
eval "$(conda shell.bash hook)"
conda activate qiime2-amplicon-2024.10

# Print commands as they run
set -x

# Trap Ctrl+C (SIGINT) and exit
trap "echo 'Script interrupted by user'; exit 1" INT


### ---- OTU level ----


# Set paths to files
FEATURE_TABLE="$INDIR/feature_table_ps99filt.biom"
METADATA="$INDIR/metadata_ps99filt.tsv"
TAXONOMY="$INDIR/taxonomy_ps99filt.tsv"



# tabulate metadata; this doesn't do anything, just allows you to visualize it in QIIME viewer
qiime metadata tabulate \
  --m-input-file $METADATA \
  --o-visualization $OUTPUT_DIR/metadata_ps99filt.qzv

# Import feature table to QIIME2
echo "Importing feature table..."
qiime tools import \
  --input-path $FEATURE_TABLE \
  --output-path $OUTPUT_DIR/feature-table_ps99filt.qza \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format

# Import taxonomy if needed
if [ -f $TAXONOMY ]; then
  echo "Importing taxonomy..."
  qiime tools import \
    --input-path $TAXONOMY \
    --output-path $OUTPUT_DIR/taxonomy_ps99filt.qza \
    --type 'FeatureData[Taxonomy]' \
    --input-format HeaderlessTSVTaxonomyFormat
fi

# collapse feature table by host species
qiime feature-table group \
  --i-table $OUTPUT_DIR/feature-table_ps99filt.qza \
  --m-metadata-file $METADATA \
  --p-axis sample \
  --m-metadata-column 'Species' \
 --p-mode mean-ceiling \
  --verbose \
  --o-grouped-table $OUTPUT_DIR/feature-table-grouped_ps99filt_meanceiling.qza

# print summary of grouped feature table for Qiime2 viewer
qiime feature-table summarize \
  --i-table $OUTPUT_DIR/feature-table-grouped_ps99filt_meanceiling.qza \
  --o-visualization $OUTPUT_DIR/feature-table-grouped_ps99filt_meanceiling.qzv

# run beta diversity on grouped feature table
qiime diversity beta-rarefaction \
  --i-table $OUTPUT_DIR/feature-table-grouped_ps99filt_meanceiling.qza \
  --p-metric 'braycurtis' \
  --p-clustering-method 'upgma' \
  --p-iterations 100 \
  --p-sampling-depth 1837 \
  --m-metadata-file $INDIR/metadata_ps99filt_grouped.tsv \
  --o-visualization $OUTPUT_DIR/beta-rarefaction_grouped_mean_1837rare.qzv

# run beta diversity on jaccard's distance (OTU level)
qiime diversity beta-rarefaction \
  --i-table $OUTPUT_DIR/feature-table-grouped_ps99filt_meanceiling.qza \
  --p-metric 'jaccard' \
  --p-clustering-method 'upgma' \
  --p-iterations 100 \
  --p-sampling-depth 1837 \
  --m-metadata-file $INDIR/metadata_ps99filt_grouped.tsv \
  --o-visualization $OUTPUT_DIR/beta-rarefaction_grouped_jaccard_mean1837rare.qzv


##### GENUS LEVEL WITH TAX GLOMMED IN R ####

# set new feature table
FEATURE_GENUS="$INDIR/feature_table_ps99filt_genus.biom"

# import feature table
qiime tools import \
  --input-path $FEATURE_GENUS \
  --output-path $OUTPUT_DIR/feature-table_ps99filt_genusinR.qza \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format

# create qzv to make sure it worked
qiime feature-table summarize \
  --i-table $OUTPUT_DIR/feature-table_ps99filt_genusinR.qza \
  --o-visualization $OUTPUT_DIR/feature-table_ps99filt_genusinR.qzv

# group by host species
qiime feature-table group \
  --i-table $OUTPUT_DIR/feature-table_ps99filt_genusinR.qza \
  --m-metadata-file $METADATA \
  --p-axis sample \
  --m-metadata-column 'Species' \
  --p-mode mean-ceiling \
  --verbose \
  --o-grouped-table $OUTPUT_DIR/feature-table-genusinR_grouped_ps99filt.qza

# create qzv to make sure it worked
qiime feature-table summarize \
  --i-table $OUTPUT_DIR/feature-table-genusinR_grouped_ps99filt.qza \
  --o-visualization $OUTPUT_DIR/feature-table-genusinR_grouped_ps99filt.qzv

# run beta diversity on bray-curtis distance (genus level)
qiime diversity beta-rarefaction \
  --i-table $OUTPUT_DIR/feature-table-genusinR_grouped_ps99filt.qza \
  --p-metric 'braycurtis' \
  --p-clustering-method 'upgma' \
  --p-iterations 100 \
  --p-sampling-depth 1275 \
  --m-metadata-file $INDIR/metadata_ps99filt_grouped.tsv \
  --o-visualization $OUTPUT_DIR/beta-rarefaction_genusinR_braycurtis_mean1275rare.qzv

# run beta diversity on jaccard's distance (genus level)
qiime diversity beta-rarefaction \
  --i-table $OUTPUT_DIR/feature-table-genusinR_grouped_ps99filt.qza \
  --p-metric 'jaccard' \
  --p-clustering-method 'upgma' \
  --p-iterations 100 \
  --p-sampling-depth 1275 \
  --m-metadata-file $INDIR/metadata_ps99filt_grouped.tsv \
  --o-visualization $OUTPUT_DIR/beta-rarefaction_genusinR_jaccard_mean1275rare.qzv
