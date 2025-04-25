### wrangling VSEARCH output to a phyloseq object
# EVS 3/2025

library(vegan)
library(phyloseq)
library(microViz)
library(tidyverse)
library(EcolUtils)
library(ggpubr)
library(ape)


### ---- function: VtoPhyloseq ----

### make function to read in VSEARCH stuff and make phyloseq
VtoPhylo <- function(path) {
  # detect OTU table
  # detect OTU table
  otabfile <- list.files(path = path, pattern = "*otutab.txt", full.names = TRUE)
  # OTU Table
  otab <- read.table(otabfile, sep = "\t", header = TRUE, comment.char = "", na.strings = "") %>% column_to_rownames(var = "X.OTU.ID")
  
  # SINTAX wrangle
  # detect SINTAX file
  sinfile <- list.files(path = path, pattern = "*sintax50.txt", full.names = TRUE)
  sin <- read.table(sinfile, sep = "\t", header = FALSE, comment.char = "", na.strings = "") %>%
    dplyr::select(V1, V4) %>%
    separate_wider_delim(cols = V4, names = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), delim = ",", too_few = "align_start") %>%
    mutate(OTU_ID = str_extract(V1, "OTU_(\\d){1,10}")) %>%
    mutate(across(!OTU_ID, ~ str_remove(.x, "[:alpha:]:"))) %>%
    column_to_rownames(var = "OTU_ID") %>%
    dplyr::select(-V1)
  
  # make phyloseq
  ps <- phyloseq(
    otu_table(otab, taxa_are_rows = TRUE),
    tax_table(sin %>% as.matrix())
  )
  
  # wrangle sample names
  sample_names(ps) <- str_remove(sample_names(ps), "^X")
  sample_names(ps) <- str_remove(sample_names(ps), "_1$")
  sample_names(ps) <- str_remove(sample_names(ps), "_S(\\d){1,4}$")
  
  # add
  sample_data(ps) <- allmeta
  
  # return phyloseq object
  return(ps)
}

## ---- make phyloseq ----

# add metadata
allmeta <- read.table("data/all_primate_metadata.txt", sep = "\t", header = TRUE) %>%
  # add the two Blanks for PRJNA686661
  full_join(data.frame(
    Run = c("Blank1_009_A09_011", "Blank2_009_A10_011"),
    Species = "BLANK",
    Accession = "PRJNA686661",
    Paper = "Sharma2022"
  )) %>%
  column_to_rownames(var = "Run")

# get phyloseq of 99% OTUs
ps99 <- VtoPhylo(path = "data") # 33K taxa
ps99 <- ps99 %>% tax_select("Fungi", "Kingdom")

### ---- rarefy and filter ----

## rarefy 99%
otab <- ps99 %>%
  otu_get() %>%
  as.data.frame()
rtab <- otab[rowSums(otab) >= 2000, ]
raretab <- rrarefy.perm(rtab, sample = 2000, n = 100, round.out = TRUE)

# re-phyloseq
rarephy99 <- phyloseq(
  otu_table(raretab, taxa_are_rows = FALSE),
  tax_table(ps99),
  sample_data(ps99)
)

# filter
filt99 <- rarephy99 %>%
  tax_filter(min_prevalence = 2, min_total_abundance = 2, min_sample_abundance = 2) 

# save filt99 object
save(filt99, file = "data/phylo_OTU99_ITS2.RData")

## for phylosymbiosis in QIIME: filter and save without rarefaction
ps99filt <- ps99 %>%
  tax_filter(min_prevalence = 2, min_total_abundance = 2, min_sample_abundance = 2)
# save
save(ps99filt, file = "data/phylo_OTU99_ITS2_filt_norarefaction.RData")

## for phylosymbiosis in QIIME at the genus level, 
# keep only OTUs with at least a family-level assignment, glom to Genus, and save
fams <- get_taxa_unique(ps99filt, "Family")
fams <- fams[!is.na(fams)]
ps99filtg <- ps99filt %>%
  tax_select(fams, "Family", strict_matches = TRUE) %>%
  tax_fix() %>%
  tax_agg("Genus")
save(ps99filtg, file = "data/phylo_OTU99_ITS2_filt_norarefaction_genus.RData")

### ---- basic stats ----

## how mnay OTUs have known genera?
k99 <- filt99 %>%
  tt_get() %>%
  as.data.frame() %>%
  drop_na(Genus) %>%
  nrow() # 2600
k99 / ntaxa(filt99) 
