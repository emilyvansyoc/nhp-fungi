## Bray-Curtis pairwise distances 
#EVS 12/2023, updated 4/2025

library(vegan)
library(phyloseq)
library(microViz)
library(tidyverse)
# to get pairwise distances into a nice dataframe:
# devtools::install_github("kylebittinger/usedist")
library(usedist)
library(cowplot)
library(rstatix)
library(ggtext) # for element markdown
library(multcompView) # for letters in boxplot
library(ggh4x) # for color strip backgrounds
library(rphylopic)
library(pairwiseAdonis)
#library(patchwork)
library(infer)
set.seed(234)


# load
load("data/phylo_OTU99_ITS2.RData")

# calculate bray
uf <- filt99 %>%
  dist_calc("bray") %>%
  microViz::dist_get()

# use 'usedist' package to convert to nice dataframe
ufdf <- usedist::dist_groups(d = uf, g = filt99@sam_data$Species)

## ---- permanovas ----

# broad anova for general summary stats
filt99 %>% dist_calc('bray') %>% dist_permanova(variables = "Species")

# pairwise anovas
sampdf <- filt99 %>% samdat_tbl()
pairwise.adonis2(uf ~ Species, data = sampdf) # sigs between all

### ---- intraspecific comparisons ----

# wrangle
dat <- ufdf %>% 
  filter(str_detect(Label, "Within")) %>% 
  mutate(newLab = Group1,
         is.within = "Within") %>% 
  dplyr::select(Item1, Item2, is.within, Distance, newLab) %>% 
  rbind(ufdf %>% 
          filter(str_detect(Label, "Between")) %>% 
          mutate(is.within = "Between") %>% 
          dplyr::select(Item1, Item2, Group1, Group2, is.within, Distance) %>% 
          pivot_longer(cols = c(Group1, Group2), names_to = "oldGroup", values_to = "newLab") %>% dplyr::select(-oldGroup)) %>%
  # make new column of names to match colors
  mutate(plotcol = case_when(
    str_detect(newLab, "cercocebus") ~ "C. agilis",
    str_detect(newLab, "pan") ~ "P. troglodytes",
    str_detect(newLab, "gorilla_gorilla") ~ "G. gorilla",
    str_detect(newLab, "gorilla_beringei") ~ "G. beringei",
    str_detect(newLab, "homo") ~ "H. sapiens",
    str_detect(newLab, "indri") ~ "I. indri",
    str_detect(newLab, "papio") ~ "P. cynocephalus",
    str_detect(newLab, "proco") ~ "P. gordonum"
  )) %>% 
  # add phylopics
  mutate(uids = case_when(
    str_detect(newLab, "indri") ~ "b28bef2c-0f7a-45e7-a815-f67cd1a98d8f",
    str_detect(newLab, "papio") ~ "72f2f854-f3cd-4666-887c-35d5c256ab0f",
    str_detect(newLab, "proco") ~ "525fde4f-596e-414f-8a2c-0837aa75792c",
    str_detect(newLab, "homo") ~ "036b96de-4bca-408e-adb6-2154fcd724ef",
    str_detect(newLab, "gorilla") ~ "142e0571-3b5f-443d-a887-b572a224ea22",
    str_detect(newLab, "pan") ~ "7133ab33-cc79-4d7c-9656-48717359abb4",
    str_detect(newLab, "cercocebus") ~ "eccbb404-c99f-41f9-8785-01a7f57f1269"
  )) %>% 
  
  # make colored panels
  mutate(fulllab = interaction(is.within, plotcol, sep = "-", lex.order = TRUE),
         is.within = factor(is.within, ordered = TRUE, levels = c("Between", "Within"))) #%>% 
# make sizes by hand
#mutate(mysize = case_when(
#Common_name == "indri" ~ 0.07,
#Common_name %in% c("olive colobus", "baboon") ~ 1,
#Common_name %in% c("gorilla beringei", "gorilla gorilla") ~ 0.07,
#Common_name %in% "human" ~ 0.05,
#Common_name %in% "agile mangabey" ~ 0.07

#))


# get increasing medians for "within" distances
dat %>%
  filter(is.within == "Within") %>%
  group_by(newLab) %>%
  get_summary_stats(Distance, type = "median_iqr") %>% 
  arrange(desc(median))

# order by decreasing medians
levs <- dat %>%
  filter(is.within == "Within") %>%
  group_by(plotcol) %>%
  summarize(med = median(Distance)) %>%
  arrange(desc(med)) %>%
  pull(plotcol)

dat$plotcol<- factor(dat$plotcol, levels = levs, ordered = TRUE)

# make dataframe of phylopics
ppic <- dat %>% 
  distinct(plotcol, uids) 

#### permutation test for intraspecific
# following: https://bookdown.org/kmbm92/Applied-Biostats/perm1.html#how-to-deal-with-nonindependence-by-permutation 

### create function to run a test for each host species
myPerm <- function(species) {
  
  # get data
  datsub<- dat %>% 
    filter(plotcol == species) 
  # calculate observed test statistic (means)
  obs_stat <- datsub %>% 
    group_by(is.within) %>% 
    summarize(mean_dis_obs = mean(Distance),
              sd_dis_obs = sd(Distance))
  diff_obs <- obs_stat %>% 
    summarize(diff = diff(mean_dis_obs)) %>% pull()
  # permute by shuffling treatment group labels to create a null distribution
  datperm <- datsub %>% 
    rep_sample_n(size = nrow(datsub), reps = 999, replace = FALSE) %>% 
    group_by(replicate) %>% 
    mutate(perm_within = sample(is.within, size = nrow(datsub), replace = FALSE)) 
  # generate the permuted test statistics
  datsum <- datperm %>% 
    group_by(replicate, perm_within) %>% 
    summarize(mean_dis = mean(Distance), .groups = "drop") 
  datdiffs <- datsum %>% 
    group_by(replicate) %>% 
    summarize(diff_dis_perm = diff(mean_dis))
  # calculate a p value for a two-sided test
  permdiffs <- datdiffs %>% 
    mutate(abs_obs_diff = abs(diff_obs),
           abs_perm_diff = abs(diff_dis_perm),
           as_or_more_extreme = abs_perm_diff >= abs_obs_diff)
  # return p value
  pval <- mean(permdiffs$as_or_more_extreme)
  return(pval)
  # plot the permuted distribution
  ggplot(permdiffs, aes(x = diff_dis_perm, fill = as_or_more_extreme)) +
    geom_histogram(color = "white") +
    scale_fill_manual(values = c("grey", "black")) +
    theme_light() %>% 
    geom_vline(xintercept = c(-1, 1)*diff_obs, color = "lightblue") +
    ggtitle(species)
  
}

#### perform test for each species
unique.species <- unique(dat$plotcol)
tests <- lapply(unique.species, myPerm)

## get summary statistics for supplementary
write.table(dat %>% 
              group_by(plotcol, is.within) %>% 
              get_summary_stats(Distance, type = "mean_ci") %>% 
              dplyr::select(plotcol, is.within, n, mean, ci) %>% 
              pivot_wider(names_from = is.within, values_from = c(n, mean, ci)),
            file = "data/supplementary_intraspecific_sumstats_meanci.txt", sep = "\t", row.names = FALSE)


### ---- interspecific comparisons: human vs others (simpler with many speices!) ----

## get data
humv <- ufdf %>% dplyr::filter(!str_detect(Label, "Within")) %>% 
  dplyr::filter(str_detect(Label, "homo_sapiens")) %>% 
  # pretty-fy
  dplyr::mutate(labplot = case_when(
    Label %in% "Between cercocebus_agilis and homo_sapiens" ~ "C. agilis",
    Label %in% "Between homo_sapiens and pan_troglodytes"  ~ "P. troglodytes",
    Label %in% "Between gorilla_gorilla and homo_sapiens" ~ "G. gorilla",
    Label %in% "Between gorilla_beringei and homo_sapiens" ~ "G. beringei",
    Label %in% "Between homo_sapiens and indri_indri" ~ "I. indri",
    Label %in% "Between homo_sapiens and papio_cynocephalus" ~ "P. cynocephalus",
    Label %in% "Between homo_sapiens and procolobus" ~ "P. gordonum"
  ))

## use letters instead of brackets
#### THIS RELIES ON ANOVA TEST FOR CODING PURPOSES BUT RESULTS ARE IDENTICAL WITH K-W AND DUNN TEST
mod <- aov(Distance ~ Label, data = humv)
t <- TukeyHSD(mod)

# save for supplementary materisl
write.table(tidy(t), file = "data/supplementary_interspecific_ANOVA.txt", sep = "\t", row.names = FALSE)

