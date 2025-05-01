# make updated colors
# EVS 4/2025

library(tidyverse)
library(RColorBrewer)
library(scales)
library(ape)

## get names
# get tree
tr <- read.tree("data/topocong_input/host_ITS2_rotated.newick")

## make plot showing clade
allmeta <- read.table("data/all_primate_metadata.txt", sep = "\t", header = TRUE)
allmeta <- allmeta %>% 
  mutate(Species = str_replace_all(Species, " ", "_")) %>% 
  mutate(Species = if_else(str_detect(Species, "troglodyte"), "pan_troglodytes", Species)) %>% 
  mutate(Species = if_else(str_detect(Species, "procolobus"), "procolobus", Species))

# add to tree
sampdf <- allmeta %>% 
  filter(Species %in% tr$tip.label) %>% 
  filter(!Accession == "PRJNA1029329") %>% 
  filter(!Accession == "PRJEB32407") %>% 
  dplyr::select(Accession, primerset, #Country, Region, Sampling_Season,
                General_phylogeny, Common_name, Family, Genus, Species) %>% 
  distinct() %>% 
  relocate(Species) %>% 
  # fix
  mutate(Common_name = case_when(
    Species %in% "gorilla_beringei"  ~ "gorilla beringei",
    Species %in% "gorilla_gorilla" ~ "gorilla gorilla",
    .default = Common_name
  )) %>% 
  mutate(plotlab = case_when(
    Species %in% "indri_indri" ~ "I. indri",
    Species %in% "papio_cynocephalus" ~ "P. cynocephalus",
    Species %in% "procolobus" ~ "P. gordonum",
    Species %in% "homo_sapiens" ~ "H. sapiens",
    Species %in% "pan_troglodytes" ~ "P. troglodytes",
    Species %in% "cercocebus_agilis" ~ "C. agilis",
    Species %in% "gorilla_beringei" ~ "G. beringei",
    Species %in% "gorilla_gorilla" ~ "G. gorilla"
  ))

# set colors
show_col(pal_brewer(palette = "Set1")(8))
ccol <- brewer_pal(palette = "Set1")(8)
#names(ccol) <- unique(sampdf$Common_name)
names(ccol) <- unique(sampdf$plotlab)
# change the yellow - too brigh
ccol[names(ccol) == "C. agilis"] <- "#FFDE21"
