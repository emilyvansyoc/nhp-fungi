## prep topo cong trees
# EVS 4/2025

library(ape)
library(tidyverse)

## ---- SET VARIABLES ----

# location to QIIME2 trees
q2 <- "data/qiime2-output/"

# location to Figree rooted trees
fg <- "data/rooted_dendrograms/"

# set output location
out <- "data/topocong_input/"

## ---- host tree ----

## HOST TREE
host <- ape::read.tree("data/timetree/my_primate_species_matched.nwk")
# keep the ITS2 tip labels
host <- keep.tip(host, tip = c("indri_indri", "gorilla_gorilla", "gorilla_beringei", "homo_sapiens", "pan_troglodytes", "papio_cynocephalus", "cercocebus_agilis", "procolobus"))

## manually root and rotate branches
h1 <- root(host, outgroup = "indri_indri")
h2 <- rotate(h1, node = 10)
h3 <- rotate(h2, node = 11)
# write out and export as PDF
write.tree(h3, file = "data/topocong_input/host_ITS2_rotated.newick")
pdf(file = "data/topocong_input/host_ITS2_rotated.pdf")
plot(h3)
dev.off()


##### ---- OTU level ----

# bray-curtis
b <- read.tree(file.path(fg, "rooted_tree_bray_OTU_meanrare1837.newick"))
#b$tip.label <- str_remove_all(b$tip.label, "'")
b1 <- rotate(b, node = 9)
b2 <- rotate(b1, node = 12)
b3 <- rotate(b2, node = 14)
#b3 <- rotate(b2, node = 13)
#b4 <- rotate(b3, node = 14)
# get node labels
bnode <- read.tree(file.path(q2, "tree_bray_OTU_meanrare1837.tre"))
b3$node.label <- bnode$node.label
write.tree(b3, file = "data/topocong_input/bray_OTU_rotated_meanrare1837.newick")
pdf(file = "data/topocong_input/bray_OTU_rotated_meanrare1837.pdf")
plot(b3)
dev.off()

# jaccard 
j <- read.tree(file.path(fg, "rooted_tree_jaccard_OTU_meanrare1837.newick"))
j1 <- rotate(j, node = 9)
j2 <- rotate(j1, node = 12)
j3 <- rotate(j2, node = 13)
j4 <- rotate(j3, node = 14)
# get node labels
jnode <- read.tree(file.path(q2, "tree_jaccard_OTU_meanrare1837.tre"))
j4$node.label <- jnode$node.label
write.tree(j4, file = file.path(out, "jaccard_OTU_rotated_meanrare1837.newick"))
pdf(file.path(out, "jaccard_OTU_rotated_meanrare1837.pdf"))
plot(j4)
dev.off()


### ---- Genus level ----

# bray
b <- read.tree(file.path(fg, "rooted_bray_genusinR_meanrare1275.newick"))
b1 <- rotate(b, node = 9)
b2 <- rotate(b1, node = 12)
# get node lables
bnode <- read.tree(file.path(q2, "tree_bray_genusinR_meanrare1275.tre"))
b2$node.label <- bnode$node.label
write.tree(b2, file.path(out, "bray_genusinR_rotated_meanrare1275.newick"))
pdf(file.path(out, "bray_genusinR_rotated_meanrare1275.pdf"))
plot(b2)
dev.off()

# jaccard
j <- read.tree(file.path(fg, "rooted_jaccard_genusinR_meanrare1275.newick"))
j1 <- rotate(j, node = 9)
j2 <- rotate(j1, node = 10)
j3 <- rotate(j2, node = 11)
j4 <- rotate(j3, node = 13)
j5 <- rotate(j4, node = 15)
# get node labels
jnode <- read.tree(file.path(q2, "tree_jaccard_genusinR_meanrare1837.tre"))
j5$node.label <- jnode$node.label
write.tree(j5, file.path(out, "jaccard_genusinR_rotated_meanrare1275.newick"))
pdf(file.path(out, "jaccard_genusinR_rotated_meanrare1275.pdf"))
plot(j5)
dev.off()
