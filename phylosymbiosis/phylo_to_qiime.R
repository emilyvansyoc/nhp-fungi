# Script to convert phyloseq object to QIIME2 format
library(phyloseq)
library(tidyverse)
library(biomformat)
library(microViz)

# define function
# Ensure row names are OTU IDs and column names are sample IDs
makeBiom <- function(ps) {
  if (taxa_are_rows(ps)) {
    # OTUs are rows, keep as is
    feature_table <- otu_table(ps)
  } else {
    # OTUs are columns, transpose
    feature_table <- t(otu_table(ps))
  }
  
  otu_table <- as.matrix(feature_table)
  otu_biom <- make_biom(data = otu_table)
  return(otu_biom)
}

## ---- write qiime2 files on the 'raw' phyloseq object without rarefaction ----


# Load the phyloseq object
load("data/phylo_OTU99_ITS2_filt_norarefaction.RData")

# remove blanks
ps99filt <- ps99filt %>% ps_filter(!Species == "BLANK")
# add underscore to species names
ps99filt <- ps99filt %>%
  ps_mutate(Species = str_replace_all(Species, " ", "_"))

# run function and save
myBiom <- makeBiom(ps99filt)
head(myBiom)
write_biom(myBiom, "data/for-qiime2/feature_table_ps99filt.biom")


# Extract sample metadata
metadata <- as.data.frame(sample_data(ps99filt))

# extract with microViz and format for qiime2
metadata <- ps99filt %>%
  samdat_tbl() %>%
  dplyr::rename(id = .sample_name)

# Write the metadata as a TSV file
write.table(metadata,
            file = "data/for-qiime2/metadata_ps99filt.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE
)

# Extract taxonomy table (optional, useful for additional analyses)
tax_table <- as.data.frame(tax_table(ps99filt))
tax_table <- tax_table %>%
  rownames_to_column(var = "Feature ID")

# Create a taxonomy string in QIIME2 format
tax_table$taxonomy <- paste(
  paste0("k__", tax_table$Kingdom),
  paste0("p__", tax_table$Phylum),
  paste0("c__", tax_table$Class),
  paste0("o__", tax_table$Order),
  paste0("f__", tax_table$Family),
  paste0("g__", tax_table$Genus),
  sep = "; "
)

# Write the taxonomy as a TSV file
write.table(tax_table[, c("Feature ID", "taxonomy")],
            file = "data/for-qiime2/taxonomy_ps99filt.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE
)

# create metadata file with only host species names
metagroup <- metadata %>%
  dplyr::select(Species) %>%
  distinct() %>%
  mutate(id = Species) %>%
  relocate(id)
# write
write.table(metagroup, "data/for-qiime2/metadata_ps99filt_grouped.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE
)

### ---- write qiime2 files on tax glommed phyloseq object ----

# load
load("data/phylo_OTU99_ITS2_filt_norarefaction_genus.RData")

myBiom <- makeBiom(ps99filtg)
head(myBiom)
write_biom(myBiom, "data/for-qiime2/feature_table_ps99filt_genus.biom")
