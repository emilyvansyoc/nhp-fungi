# Sequence divergence analysis and figure - final
# EVA 11/2024

library(multcompView)
library(ape)
library(tidyverse)
library(usedist)
library(ggpubr)
library(ggtext)
library(cowplot)
library(vegan)
#library(phyloseq)
#library(microViz)
library(cowplot)

# ---- make function ----


myDiv <- function(alignment, # fasta file of genus level alignment
                  OTU, # character of OTU to pull
                  dist.type # from dist.dna 
) {
  
  # print progress message
  cat("\n ##### Working on sequence divergence for ", OTU, " ######## \n")
  
  # subset alignment to get the OTU of interest
  sub <- alignment[str_detect(labels(alignment), paste0(OTU, ";")), ]
  # print synopsis to screen
  print(sub)
  
  # get the divergence
  div <- dist.dna(sub, model = dist.type)
  
  # data frame-ize
  groupdf <- data.frame(ids = labels(div)) %>% 
    mutate(hominid = case_when(
      str_detect(ids, "MountGorilla") ~ "G. beringei",
      str_detect(ids, "LowGorilla") ~ "G. gorilla",
      str_detect(ids, "Human") ~ "H. sapiens",
      str_detect(ids, "Chimp") ~ "P. troglodytes"
    ))
  disdf <- dist_groups(div, g = groupdf$hominid)
  
  # print the number of exact sequence variants for each hominid
  cat("\n Number of sequence variants for each hominid: \n")
  print(groupdf %>% group_by(hominid) %>% count())
  
  # set up for plot
  forplot <- disdf %>% 
    filter(!str_detect(Label, "Within")) %>% 
    mutate(labplot = case_when(
      Label == "Between G. gorilla and P. troglodytes" ~ "Between *P. troglodytes* and *G. gorilla*",
      Label == "Between H. sapiens and P. troglodytes" ~ "Between *H. sapiens* and *P. troglodytes*",
      Label == "Between G. beringei and P. troglodytes" ~ "Between *P. troglodytes* and *G. beringei*",
      Label == "Between G. gorilla and H. sapiens" ~ "Between *H. sapiens* and *G. gorilla*",
      Label == "Between G. beringei and G. gorilla" ~ "Between *G. gorilla* and *G. beringei*",
      Label == "Between G. beringei and H. sapiens" ~ "Between *H. sapiens* and *G. beringei*"
    )) 
  
  # make test
  mod <- aov(Distance ~ labplot, data = forplot)
  t <- TukeyHSD(mod)
  cld <- multcompLetters4(mod, t, reversed = TRUE)
  cldf <- as.data.frame.list(cld$labplot) %>% 
    rownames_to_column(var = "labplot") %>% 
    dplyr::select(labplot, Letters) %>% 
    # get positions
    full_join(forplot %>% group_by(labplot) %>% get_summary_stats(Distance, type = "common") %>% dplyr::select(mean, sd, labplot)) 
  # print stat test
  cat("\n ANOVA test: \n")
  print(cldf)
  
  # plot
  #myplot <- ggplot(data = forplot, 
  #                aes(x = fct_reorder(labplot, desc(Distance)), y = Distance)) +
  
  # geom_point(position = position_jitter(width = 0.2, seed = 123), alpha = 0.5, color = "darkgrey") +
  # geom_boxplot(size = 1.5, outlier.shape = NA, notch = FALSE, fill = "lightcyan4", alpha = 0.7) +
  # geom_text(data = cldf, aes(x = labplot, label = Letters, y = 0),
  # size = 10) +
  #ylim(-0.1, 1.1) +
  #scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  #coord_flip() +
  # theme_pubr() +
  #labs(x = "", y = "Nucleotide sequence divergence", title = OTU) +
  #theme(axis.text.y = element_markdown())
  
  # plot in the window
  #plot(myplot)
  
  # return stat test and plot
  #return(list(forplot, myplot))
  return(list(forplot))
  
  
  
}

### ---- run on cophylogenetic fungi ----

# get parafit/PACo results for significant fungi
load("cophylogeny/cophylogeny_scripts/data/all_genera_list.RData")

# get all genera and OTUs
allgen <- read.table("cophylogeny/cophylogeny_scripts/data/all_pacoparafit_results_wr2.txt", sep = "\t", header = TRUE)

## ---- divergence: raw nucleotide differences ----


# loop through and get divergence estimates in a big dataframe
myotus <- unique(allgen$otu)
alldf <- data.frame()
for(i in 1:length(myotus)) {
  
  # get genus name
  sub <- allgen %>% filter(otu == myotus[i])
  genname <- unique(sub$labname)
  # loop
  out <- myDiv(alignment = read.dna(paste0("/Users/epb5360/bioinformatics/r-projects/hominid/updated_methods/wmangabey/codiv_wmangabey/align_genus_renamed//unique_clipkit_derep_haplotypes_", genname, ".fasta"), format = "fasta"),
               OTU = myotus[i],
               dist.type = "N")
  
  # add to dataframe
  alldf <- rbind(alldf, out[[1]] %>% mutate(otu = myotus[i], Genus = genname))
}

# get summary stats for total sequence divergence across all OTUs
alldf %>% group_by(otu, Genus) %>% get_summary_stats(Distance, type = "common") %>% arrange(desc(mean)) %>% print(n=45)

## ---- divergence: proportions ----


propdf <- data.frame()
for(i in 1:length(myotus)) {
  
  # get genus name
  sub <- allgen %>% filter(otu == myotus[i])
  genname <- unique(sub$labname)
  # loop
  out <- myDiv(alignment = read.dna(paste0("/Users/epb5360/bioinformatics/r-projects/hominid/updated_methods/wmangabey/codiv_wmangabey/align_genus_renamed//unique_clipkit_derep_haplotypes_", genname, ".fasta"), format = "fasta"),
               OTU = myotus[i],
               dist.type = "raw")
  
  # add to dataframe
  propdf <- rbind(propdf, out[[1]] %>% mutate(otu = myotus[i], Genus = genname))
}

## ---- use proportional divergence to calibrate clock ----


# start with averages and use confidence intervals 
avgdf2 <- propdf %>% 
  group_by(Label, otu, Genus) %>% 
  get_summary_stats(Distance, type = "mean_sd") %>% 
  left_join(allgen, by = "otu") %>% drop_na(pacor2)

# using the average between HS and PT to 'calibrate' for each OTU
hs.pt1 <- avgdf2 %>% 
  ungroup() %>% 
  filter(Label %in% c("Between H. sapiens and P. troglodytes")) %>% 
  # assume that host speciation is 6MYA
  mutate(hspt.div = 6) %>% 
  mutate(hspt.clock = mean / hspt.div) %>% # mutation rate every 1MY 
  dplyr::select(otu, Genus, hspt.clock)
summary(hs.pt1$hspt.clock) # mean 0.2% every MYA

# add this back to the dataset and use it to calculate the other 'clocks'
clockdf2 <- propdf %>% 
  ungroup() %>% 
  left_join(allgen, by = "otu") %>% drop_na(pacor2) %>% 
  filter(issig == TRUE) %>% 
  filter(!Label %in% c("Between H. sapiens and P. troglodytes")) %>% 
  left_join(hs.pt1) %>% 
  # calculate the divergence time of other speciations
  mutate(divtime = Distance / hspt.clock)

# get confidence intervals around the estimates
clockest <- clockdf2 %>% 
  group_by(otu, Genus, Species, Label) %>% 
  get_summary_stats(divtime, type = "full") %>% 
  mutate(cihi = mean + ci,
         cilo = mean - ci)

### wrangle and save for summary statistics
ssum <- clockest %>% 
  # make column of pretty-fied names that match codiv fig
  mutate(taxname = str_replace(Species, " Genus", " spp."),
         taxname = str_replace(taxname, " Order", " spp."),
         taxname = str_replace(taxname, "_", " ")) %>% 
  mutate(taxname = if_else( otu == "OTU_456", "Xylaria spp. (1)", taxname),
         taxname = if_else(otu == "OTU_384", "Xylaria spp. (2)", taxname))  %>% 
  dplyr::select(Label, otu, taxname, mean, ci)


### ---- save/load image ----

# to execute parts of this script, load the image 
#save.image(file = "data/image_sequencedivergence.RData")

