### parafit and PACO results
# EVS 4/2024, updated 5/2/2024 with latest PACo OTU results from ROAR

library(ape)
library(vegan)
library(tidyverse)
library(phyloseq)
library(microViz)
library(Biostrings)
library(purrr)

## ---- OTU ----

## get trees
load("cophylogeny/cophylogeny_scripts/data/all_otu_subtrees.RData")

# get Parafit results
load( "cophylogeny/cophylogeny_scripts/data/parafit_results_OTU.RData")
# get sigs
paradf <- paradf %>% 
  mutate(pfdr = p.adjust(pval, method = "fdr"),
         pbon = p.adjust(pval, method = "bonferroni"))

## get PACo output: 
dir <- "cophylogeny/cophylogeny_scripts/data/PACo_OTU_quasi_p100_10-15-2024/"
myfi <- list.files(dir, pattern = ".txt", full.names = TRUE)
#myfi <- myfi[!str_detect(myfi, "p10_")]
pacodf <- data.frame()
for(i in 1:length(myfi)) {
  fi <- read.table(myfi[i], sep = "\t", header = TRUE)
  pacodf <- rbind(pacodf, fi)
}
pacodf <- pacodf %>% 
  mutate(pfdr = p.adjust(paco.p, method = "fdr"),
         pbon = p.adjust(paco.p, method = "bonferroni"))

## join
allotu <- paradf %>%
  dplyr::select(otu, para.p = pval, parastat, para.pfdr = pfdr, para.pbon = pbon) %>% 
  full_join(pacodf %>% dplyr::select(otu = subtree, paco.ss, paco.p, paco.pfdr = pfdr, paco.pbon = pbon)) %>% 
  # calculate paco R2
  mutate(pacor2 = 1-paco.ss)
summary(allotu$pacor2)
## get sigs from both
bsigs <- allotu %>% filter(para.pfdr < 0.05 & paco.pfdr < 0.05)
# for now; 11
summary(bsigs$pacor2)
summary(bsigs$parastat) # not as easily interpretable as the r2
## 
# get biggest effect size
bsigs %>% filter(pacor2 > 0.25) # only OTU_144 
quantile(bsigs$pacor2)

# save for supplementary materials
# get phyloseq
load("data/cophylogeny_data/hominidphylo_ITS.RData")

tax <- psf %>% tax_fix() %>% tax_select(allotu$otu, strict_matches = TRUE, n_typos = 0) %>% 
  tt_get() %>% as.data.frame() %>% 
  rownames_to_column(var = "otu") %>% 
  full_join(allotu) %>% 
  dplyr::select(otu, Phylum, Class, Order, Family, Genus, Species,
                parastat, para.pfdr, pacor2, paco.pfdr) %>%
  arrange(desc(pacor2))



### ----- get fungal information -----



## get the rest of the OTUs in bsigs
sigotus <- bsigs$otu
sigps <- psf %>% tax_fix() %>% tax_select(sigotus, strict_matches = TRUE, n_typos = 0) %>% 
  tt_get() %>% as.data.frame() %>% 
  rownames_to_column(var = "otu") %>% full_join(bsigs) %>% arrange(desc(pacor2))


### ---- what is the prevalence of significant otus ----

pa <- psf %>% ps_filter(!Group == "Mangabey") %>% 
  
  ps_mutate(SciName = case_when(
    SpeciesCaptive %in% "Wild_Chimp" ~ "P_troglodytes_schweinfurthii",
    SpeciesCaptive %in% "Wild_Lowland Gorilla" ~ "Gorilla_gorilla",
    SpeciesCaptive %in% "Wild_Mountain Gorilla" ~ "Gorilla_beringei",
    SpeciesCaptive %in% "Human_Bantu" ~ "Homo_sapien",
    SpeciesCaptive %in% "Human_BaAka" ~ "Homo_sapien"
  )) %>% 
  tax_select(sigps$otu) %>% tax_transform('pa') %>% ps_melt() %>% 
  group_by(OTU, SciName) %>% summarize(tot = sum(Abundance))

# get species totals to calculate prevalence
stot <- psf %>% ps_filter(!Group == "Mangabey") %>% 
  
  ps_mutate(SciName = case_when(
    SpeciesCaptive %in% "Wild_Chimp" ~ "P_troglodytes_schweinfurthii",
    SpeciesCaptive %in% "Wild_Lowland Gorilla" ~ "Gorilla_gorilla",
    SpeciesCaptive %in% "Wild_Mountain Gorilla" ~ "Gorilla_beringei",
    SpeciesCaptive %in% "Human_Bantu" ~ "Homo_sapien",
    SpeciesCaptive %in% "Human_BaAka" ~ "Homo_sapien"
  )) %>% samdat_tbl() %>% group_by(SciName) %>% count()
# calculate prevalence
pa <- pa %>% mutate(prevalence = case_when(
  SciName == "Gorilla_beringei" ~ tot / stot$n[stot$SciName == "Gorilla_beringei"],
  SciName == "Gorilla_gorilla" ~ tot / stot$n[stot$SciName == "Gorilla_gorilla"],
  SciName == "P_troglodytes_schweinfurthii" ~ tot / stot$n[stot$SciName == "P_troglodytes_schweinfurthii"],
  SciName == "Homo_sapien" ~ tot / stot$n[stot$SciName == "Homo_sapien"]
))

### calculate the relative abundance of cophylo OTUs
rel <- psf %>% tax_transform("compositional") %>% tax_select(sigotus, strict_matches= TRUE) %>% 
  ps_melt()

# get the summary statistics of these
rel %>% group_by(OTU, Genus, Species) %>% get_summary_stats(Abundance, type = "mean_sd")
