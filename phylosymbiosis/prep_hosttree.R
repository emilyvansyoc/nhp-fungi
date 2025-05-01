# hominid phylogeny from timetree
# EVS 3/2025

library(ape)
library(tidyverse)

# read in timetree
tree <- read.tree("data/timetree/primate_species.nwk")
# 449 tips; rooted

# get tip labels
labs <- str_to_lower(tree$tip.label)
labs <- str_replace_all(labs, "_", " ")

# read in metadata
meta <- read.table("data/all_primate_metadata.txt", header = TRUE, sep = "\t")

# get host species
sp <- unique(meta$Species)

# make exact matches
match <- labs[labs %in% sp]

# get non-matches
nomatch <- sp[!sp %in% match]
# pan troglodytes schweinfurthii is OK - subspecies
# the only procolobus in the tree is "procolobus verus"
# the only phaner in the tree is just "phaner"
# maybe - OK to use the two genera for these?
meta %>%
  filter(Genus == "procolobus" | Genus == "phaner") %>%
  select(Species) %>%
  distinct()
# both are the only species representation in the genus

# force the matches
meta1 <- meta %>%
  mutate(Species = case_when(
    Species %in% "pan troglodytes schweinfurthii" ~ "pan troglodytes",
    Species %in% "procolobus gordonorum" ~ "procolobus",
    Species %in% "phaner pallescens" ~ "phaner",
    TRUE ~ Species
  ))
# re-match
match1 <- labs[labs %in% meta1$Species]
# get the non-matches
nomatch1 <- unique(meta1$Species)[!unique(meta1$Species) %in% match1]
# fix this one by hand

tree1 <- tree
tree1$tip.label <- str_replace_all(str_to_lower(tree1$tip.label), "_", " ")
tree1$tip.label <- if_else(str_detect(tree1$tip.label, "procolobus"), "procolobus", tree1$tip.label)

# prune
tree2 <- keep.tip(tree1, unique(meta1$Species)) # 33 - looks good

plot(tree2)
tiplabels()

# save
write.tree(tree2, file = "data/timetree/my_primate_species_matched.nwk")
