## prep topo cong trees
# EVS 4/2025

library(ape)
library(tidyverse)

## ---- SET VARIABLES ----

# location to QIIME2 trees
q2 <- "all_primate/qiime2/output/output_trees/"

# location to Figree rooted trees
fg <- "all_primate/qiime2/output/output_trees/FigTree_output/"

# set output location
out <- "all_primate/qiime2_topocong/input/"

## ---- host tree ----

## HOST TREE
host <- ape::read.tree("all_primate/primate_tree/timetree/my_primate_species_matched.nwk")
# keep the ITS2 tip labels
host <- keep.tip(host, tip = c("indri_indri", "gorilla_gorilla", "gorilla_beringei", "homo_sapiens", "pan_troglodytes", "papio_cynocephalus", "cercocebus_agilis", "procolobus"))

## manually root and rotate branches
h1 <- root(host, outgroup = "indri_indri")
h2 <- rotate(h1, node = 10)
h3 <- rotate(h2, node = 11)
# write out and export as PDF
write.tree(h3, file = "all_primate/qiime2_topocong/input/host_ITS2_rotated.newick")
pdf(file = "all_primate/qiime2_topocong/input/host_ITS2_rotated.pdf")
plot(h3)
dev.off()

### get host tree for just hominid
h4 <- keep.tip(h3, tip = c("gorilla_gorilla", "gorilla_beringei", "homo_sapiens", "pan_troglodytes", "cercocebus_agilis"))
h4 <- root(h4, outgroup = "cercocebus_agilis")
write.tree(h4, file.path(out, "host_ITS2_hominid_rotated.newick"))
pdf(file.path(out, "host_ITS_hominid_rotated.pdf"))
plot(h4)
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
write.tree(b3, file = "all_primate/qiime2_topocong/input/bray_OTU_rotated_meanrare1837.newick")
pdf(file = "all_primate/qiime2_topocong/input/bray_OTU_rotated_meanrare1837.pdf")
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

# aitchison
a <- read.tree(file.path(fg, "rooted_tree_aitchison_OTU_meanrare1837.newick"))
a1 <- rotate(a, node = 9)
a2 <- rotate(a1, node = 12)
a3 <- rotate(a2, node = 14)
# get node labels
anode <- read.tree(file.path(q2, "tree_aitchison_OTU_meanrare1837.tre"))
a3$node.label <- anode$node.label
write.tree(a3, file.path(out, "aitchison_OTU_rotated_meanrare1837.newick"))
pdf(file.path(out, "aitchison_OTU_rotated_meanrare1837.pdf"))
plot(a3)
dev.off()

# canberra
c <- read.tree(file.path(fg, "rooted_tree_canberra_OTU_meanrare1837.newick"))
c1 <- rotate(c, node = 9)
c2 <- rotate(c1, node = 11)
c3 <- rotate(c2, node = 12)
c4 <- rotate(c3, node = 13)
# get node labels
cnode <- read.tree(file.path(q2, "tree_canberra_OTU_meanrare1837.tre"))
c4$node.label <- cnode$node.label
write.tree(c4, file.path(out, "canberra_OTU_rotated_meanrare1837.newick"))
pdf(file.path(out, "canberra_OTU_rotated_mean1837.pdf"))
plot(c4)
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

# aitchison
a <- read.tree(file.path(fg, "rooted_aitchison_genusinR_meanrare1275.newick"))
a1 <- rotate(a, node = 9)
# get node labels
anode <- read.tree(file.path(q2, "tree_aitchison_genusinR_meanrare1275.tre"))
a1$node.label <- anode$node.label
write.tree(a1, file.path(out, "aitchison_genusinR_rotated_meanrare1275.newick"))
pdf(file.path(out, "aitchison_genusinR_rotated_meanrare1275.pdf"))
plot(a1)
dev.off()

# canberra
c <- read.tree(file.path(fg, "rooted_canberra_genusinR_meanrare1275.newick"))
c1 <- rotate(c, node = 9)
c2 <- rotate(c1, node = 11)
c3 <- rotate(c2, node = 14)
c4 <- rotate(c3, node = 12)
# get node labels
cnode <- read.tree(file.path(q2, "tree_canberra_genusinR_meanrare1275.tre"))
c4$node.label <- cnode$node.label
write.tree(c4, file.path(out, "canberra_genusinR_rotated_meanrare1275.newick"))
pdf(file.path(out, "canberra_genusinR_rotated_meanrare1275.pdf"))
plot(c4)
dev.off()

## ---- OTU level hominid ----

# read tree - bray
b <- read.tree(file.path(fg, "rooted_bray_OTU_hominid_meanrare5733.newick"))
b1 <- rotate(b, node = 6)
# get node labels
bnode <- read.tree(file.path(q2, "tree_bray_OTU_hominid_meanrare5733.tre"))
b1$node.label <- bnode$node.label
write.tree(b1, file.path(out, "bray_OTU_hominid_rotated_meanrare5733.newick"))
pdf(file.path(out, "bray_OTU_hominid_rotated_meanrare5733.newick"))
plot(b1)
dev.off()

## ---- OTU level with merged reads ----

# read three - bray
b <- read.tree(file.path(fg, "rooted_bray_OTU_merged_meanrare1476.newick"))
b1 <- rotate(b, node = 9)
b2 <- rotate(b1, node = 13)
b3 <- rotate(b2, node = 15)
# get node labels
bnode <- read.tree(file.path(q2, "tree_bray_OTU_merged_meanrare1476.tre"))
b3$node.label <- bnode$node.label
write.tree(b3, file.path(out, "bray_OTU_merged_rotated_meanrare1476.newick"))
pdf(file.path(out, "bray_OTU_merged_rotated_meanrare1476.newick"))
plot(b3)
dev.off()
