#### run cophylogeny scripts

library(Biostrings)
library(tidyverse)
library(phyloseq)
library(microViz)
library(ggpubr)
library(ape)
library(tidyverse)
library(doParallel)
library(foreach)
library(parallel)

## load phyloseq
load("data/cophylogeny_data/hominidphylo_ITS.RData")
# remove Mangabey 
psfix <- psf %>% tax_fix() %>% 
  ps_filter(!Group == "Mangabey")

# get genera
gen <- get_taxa_unique(psfix, "Genus")

# get names of fungal species
tt <- psfix %>% tt_get() %>% as.data.frame() %>% 
  mutate(Species = str_replace_all(Species, " ", "-")) %>% 
  rownames_to_column(var = "OTUID")

## get file to match each read to its assigned OTU
# this is output from VSEARCH
uc <- read.table("data/cophylogeny_data/usearch.uc", sep = "\t", header = FALSE) %>% 
  dplyr::select(V1, V9, V10)
names(uc) <- c("HIT", "READ", "OTU")
dat <- uc %>% 
  filter(HIT == "H") %>% 
  mutate(otulab = str_extract(OTU, "OTU_(\\d){1,10}"))

# get the file of all sequences
allfa <- readDNAStringSet("data/cophylogeny_data/all.fasta") ### not uploaded to Github (a little bit big)

# get metadata for hominid name
# this was shared by the study authors
meta <- readxl::read_xlsx("restricted/EVS_MetadataITS.xlsx") %>% 
  mutate(hominid = case_when(
    str_detect(Group, "Human") ~ "Human",
    str_detect(Group, "Chimp") ~ "Chimp",
    str_detect(Group, "Lowland Gorilla") ~ "LowGorilla",
    str_detect(Group, "Mountain Gorilla") ~ "MountGorilla"
  ))

# create a "master list" of sequence IDs, OTUs, and hominid species
ids <- names(allfa)
iddf <- data.frame(fullname = ids) %>% 
  mutate(FileID = str_remove(fullname, "_S(\\d){1,5}\\.(\\d){1,5};size=(\\d){1,10}")) %>% 
  # join metadata to get hominid
  left_join(meta %>% dplyr::select(FileID, hominid)) %>% 
  # remove NA seqs (Mangabey)
  drop_na(hominid) %>% 
  # join seqs to get OTU
  left_join(dat, by = c("fullname" = "READ")) %>% 
  # join tax table to get fungi name
  left_join(tt %>% dplyr::select(OTUID, Species), by = c("otulab" = "OTUID")) %>% 
  # drop NA OTUs (some OTUs were removed during filtering etc)
  drop_na(Species) %>% 
  # make new names for output files
  mutate(halfname = str_remove(fullname, ";size=(\\d){1,5}")) %>% 
  mutate(newids = paste(halfname, otulab, Species, hominid, sep = ";")) %>% 
  mutate(halfid = paste(otulab, Species, hominid, sep = ";"))


# ---- make function to get fungal exact sequence variants ('haplotypes') ----

myHaplo <- function(genus) {
  
  # subset the phyloseq
  sub <- psfix %>% tax_select(genus, "Genus", strict_matches = TRUE, n_typos = 0)
  
  # get OTUs
  otus <- taxa_names(sub)
  
  # subset list of IDs
  newids <- iddf %>% filter(otulab %in% otus)
  
  ### similarly to OTU-level, ensure each hominid is represented at least once (otherwise can't have codiv and are wasting comparisons!) 
  if(length(unique(newids$hominid)) < 4) {
    message("\n Genus ", genus, " does not have all hominids; it has only ", unique(newids$hominid), " ...skipping... \n")
  } else {
    
    # if all four hominids are represented, proceed
    
    # subset sequences
    seqs <- allfa[names(allfa) %in% newids$fullname ]
    # sort to match (ugh)
    seqs <- seqs[match(names(seqs), newids$fullname)]
    
    # rename sequences with fungi, OTU name, and hominid
    names(seqs) <- newids$halfid
    
    # remove spaces from family
    famname <- gsub(" ", "_", genus)
    
    # write out to file
    writeXStringSet(seqs, filepath = paste0(outdir, "/haplotypes_", famname, ".fasta"),
                    append = FALSE)
    
  }
  
}

### ---- run for haplotypes ----

# set output directory
outdir <- "cophylogeny/data/genus_haplotypes/" #### not run directly in Github; this creates a bunch of files some of which are a little large 
if(!dir.exists(outdir)) {dir.create(outdir)}

# run function
lapply(gen, function(x) myHaplo(x))

## ---- run other functions through command line ----

setwd("cophylogeny/cophylogeny_scripts/")

# dereplicate
system("bash run-vsearch-derep.sh data/genus_haplotypes data/genus_derep/")

# align
system("bash mafft-alignments.sh data/genus_derep/ data/align_genus/")

# rename
system("bash run-seqkit-rename.sh data/align_genus/ data/genus_renamed/")

## ----- make NJ trees ----

### make neighbor-joining and UPGMA trees for alignments
# make genus-level alignments


dir <- "data/genus_renamed/"
myfi <- list.files(dir, pattern = "clipkit", full.names = TRUE)

## remove genera that cannot be aligned with ITS: Fusarium, Asperillus, Penicillum, Cortinarius

# Cortinarius are not in the dataset
myfi <- myfi[!str_detect(myfi, "unique_clipkit_derep_haplotypes_Aspergillus.fasta")]
myfi <- myfi[!str_detect(myfi, "unique_clipkit_derep_haplotypes_Fusarium.fasta")]
myfi <- myfi[!str_detect(myfi, "unique_clipkit_derep_haplotypes_Penicillium.fasta")]


## output directories to write tree files
njdir <- "data/genus_trees_NJ/"
if(!dir.exists(njdir)) {dir.create(njdir)}


### RUN LOOP IN PARALLEL
n.cores <- parallel::detectCores() - 2
# create a backend
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK",
  outfile = "tmp/make-genus-trees-stdout.log" # prints stdout
)

#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
#check if it is registered (optional)
foreach::getDoParRegistered()

#### RUN LOOP
foreach(i = 1:length(myfi),
        .errorhandling = "pass" # I have my own error catch below so errors in the loop can be passed 
) %dopar% {
  
  myname <- str_remove(myfi[i], "data/unique_clipkit_derep_haplotypes_")
  myname <- str_remove(myname, ".fasta")
  
  # print progress
  cat(" \n working on file ", myname, " ... \n")
  
  # read in alignment
  aln <- read.dna(myfi[i], format = "fasta")
  
  # create distance matrix with Jukes Cantor model
  # pairwise deletion removes the site with missing data for all sequences (true gaps)
  # this is necessary because tree building will fail
  dist <- dist.dna(aln, model = "JC69", pairwise.deletion = TRUE)
  
  tryCatch(
    
    {
      # create neighbor joining tree
      # use "njs" to allow missing values 
      tree <- njs(dist)
      
      # write to file
      ## pretty-fy file name
      
      write.tree(tree, file = paste0(njdir, "/NJtree_", myname, ".newick"))
      #write.tree(dend, file = paste0(dgdir, "/UPGMAtree_", myname, ".newick"))
      cat("\n done with tree ", myname, "\n")
    },
    error = function(e) {cat("\n Alignment", myfi[i], " cannot form one or both trees \n")}
    
  ) # close tryCatch
  
  
} # close parallel loop

### close cluster
parallel::stopCluster(cl = my.cluster)

## ---- fix negative branch lengths ----

## output directory to 'fixed' trees
outdir <- "data/genus_trees_NJ_nonegbl/"
if(!dir.exists(outdir)) {dir.create(outdir)}

## get trees
myfi <- list.files(njdir, full.names = TRUE)

## loop through
for(i in 1:length(myfi)) {
  
  myname <- str_remove(myfi[i], "data/NJtree_")
  myname <- str_remove(myname, ".newick")
  
  # read tree
  tr <- read.tree(myfi[i])
  
  # replace negative branch edges
  tr$edge.length[tr$edge.length < 0] <- 0
  
  # resolve these now zero lengths into polytomies
  td <- di2multi(tr)
  
  # print how many internal nodes were resolved
  cat("\n tree", myname, "resolved", tr$Nnode - td$Nnode, "internal nodes \n ")
  
  # write out the resolved tree
  write.tree(td, file = paste0(outdir, "resolved_", myname, ".newick"))
  
  
}

## ---- get OTU subtrees ----

## directory to trees
tdir <-  "data/genus_trees_NJ_nonegbl/"
gtrees <- list.files(tdir, pattern = ".newick", full.names = TRUE)
allgtrees <- list()

# read each tree
alltrees <- list()
for(i in 1:length(gtrees)) {
  
  # get tree
  tr <- read.tree(gtrees[i])
  allgtrees[[i]] <- tr
  
  # get tip labels
  tips <- tr$tip.label
  
  # get OTUs
  otus <- data.frame(
    tiplabs = tr$tip.label
  ) %>% 
    mutate(OTU = str_extract(tiplabs, "OTU_(\\d){1,10}"),
           tip1 = str_remove(tiplabs, "_(\\d){1,10}$"),
           hominid = str_extract(tip1, "[:alpha:]+$"))
  num <- unique(otus$OTU)
  
  subtrees <- list()
  for(j in 1:length(num)) {
    
    # subset
    sub <- otus %>% filter(OTU == num[j])
    # get all 4 hominids
    if(length(unique(sub$hominid)) == 4) {
      
      # subset tree
      subtr <- keep.tip(tr, tip = sub$tiplabs)
      
      # add to list for paco
      subtrees[[j]] <- subtr
      
      names(subtrees)[[j]] <- num[j]
    } # close if statement 
    #else {subtrees[[j]] <- NULL}
    
  } # close j loop
  
  # remove NULL trees
  subtrees <- subtrees[lengths(subtrees) > 0]
  
  ## combine "big" list
  alltrees <- c(alltrees, subtrees)
  
} # close i loop


## save for PACo run
save(alltrees, file = "data/all_otu_subtrees.RData")

## ---- run parafit ----

# runs parafit on all OTU trees
source("run-parafit.R")

## ---- run PacCo ----

### this is set up to run as a job submission on a SLURM cluster
# low resource usage but the bigger trees run for weeks 
# source("run-PacCo.R")



