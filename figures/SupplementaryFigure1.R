### topo cong plots for supplementary materials
# EVS 4/2025


library(ggpubr)
library(ape)
library(tidyverse)
library(ggtree)
library(treeio)
library(ggtext)
library(gridExtra)
library(rphylopic)
library(ggimage)

# set colors
library(RColorBrewer)
library(scales)


# get phylopics- old images
h.img <- get_phylopic(uuid = "036b96de-4bca-408e-adb6-2154fcd724ef", preview = TRUE)
g.img <- get_phylopic(uuid = "142e0571-3b5f-443d-a887-b572a224ea22", preview = TRUE)
c.img <- get_phylopic(uuid = "7133ab33-cc79-4d7c-9656-48717359abb4", preview = TRUE)
# 5/8/2024; get phylopic for fungi
f.img <- get_phylopic(uuid = "aaecd181-feb8-4203-8c64-f46384257e59", preview = TRUE) # hyphae
f.img1 <- get_phylopic(uuid = "e602729e-044f-4d2b-bab2-64a87a0b48c7", preview = TRUE) # yeast

# add new images
pap.img <- get_phylopic(uuid = "72f2f854-f3cd-4666-887c-35d5c256ab0f", preview = TRUE)
proc.img <- get_phylopic(uuid = "525fde4f-596e-414f-8a2c-0837aa75792c", preview = TRUE)
ind.img <- get_phylopic(uuid = "b28bef2c-0f7a-45e7-a815-f67cd1a98d8f", preview = TRUE)
m.img <- get_phylopic(uuid = "eccbb404-c99f-41f9-8785-01a7f57f1269", preview = TRUE)

### ---- wrangle ----

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
  dplyr::select(Accession, primerset,# Country, Region, Sampling_Season,
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


### ---- make host tree phylopic ----

# convert to tree data object
host.td <- tr %>% as_tibble() %>% as.treedata()
# wrangle
host.td <- host.td %>% 
  left_join(sampdf %>% dplyr::rename(label = Species)) %>% 
  # add phylopics
  mutate(uids = case_when(
    Common_name %in% "indri" ~ "b28bef2c-0f7a-45e7-a815-f67cd1a98d8f",
    Common_name %in% "baboon" ~ "72f2f854-f3cd-4666-887c-35d5c256ab0f",
    str_detect(Common_name, "colobus") ~ "525fde4f-596e-414f-8a2c-0837aa75792c",
    str_detect(Common_name, "human") ~ "036b96de-4bca-408e-adb6-2154fcd724ef",
    str_detect(Common_name, "gorilla") ~ "142e0571-3b5f-443d-a887-b572a224ea22",
    str_detect(Common_name, "chimp") ~ "7133ab33-cc79-4d7c-9656-48717359abb4",
    str_detect(Common_name, "mangabey") ~ "eccbb404-c99f-41f9-8785-01a7f57f1269"
  )) %>% 
  # make sizes by hand
  mutate(mysize = case_when(
    Common_name == "indri" ~ 0.07,
    Common_name %in% c("olive colobus", "baboon") ~ 1,
    Common_name %in% c("gorilla beringei", "gorilla gorilla") ~ 0.07,
    Common_name %in% "human" ~ 0.05,
    Common_name %in% "agile mangabey" ~ 0.07
    
  ))

# sizes
csize <- c(0.07, 0.1, 0.08, 0.05, 0.07, 0.09, 0.07, 0.07)
names(csize) <- unique(sampdf$plotlab)

## make plot
hostfig <- ggtree(host.td, size = 1.5, ladderize = FALSE) +
  geom_rootedge(rootedge = 5, size = 1.5) +
  #geom_hilight(node = 13, type = "gradient", gradient.direction = "tr") +
  geom_tiplab(aes(image = uids, color = plotlab, size = plotlab), offset = 2, geom = "phylopic") +
  #geom_tippoint(aes(color = plotlab, shape = Family), size = 7) +
  #theme_bw() + 
  geom_treescale(x = 0, y = 0, offset = 0.1) + 
  xlim(-10, 80) +
  theme(#legend.text = element_text(size = 16),
    #legend.title = element_text(size = 20),
    legend.position = "top",
    plot.margin = margin(l = 0,r=0),
    plot.background = element_rect(fill = "white", color = NULL),
    panel.background = element_rect(fill = "white", color = NULL),
    legend.background = element_rect(fill = "white", color = NULL),
    legend.text = element_text(face = "italic", size = 12),
    legend.title = element_blank()) +
  scale_color_manual(values = ccol) +
  scale_size_manual(values = csize)

hostfig.nolegend <- hostfig + theme(legend.position = "none")


### ---- fungal dendrogram: OTU level, Jaccard ----

# read in
fd <- read.tree("data/topocong_input/jaccard_OTU_rotated_meanrare1837.newick")
# fix tip labels
fd$tip.label <- if_else(str_detect(fd$tip.label, "troglody"), "pan_troglodytes", fd$tip.label)
fd$tip.label <- if_else(str_detect(fd$tip.label, "procolobus"), "procolobus", fd$tip.label)

# convert to treedata
treetib <- as_tibble(fd) %>% as.treedata()

# wrangle
treetib <- treetib %>% 
  mutate(nodelabs = if_else(isTip == TRUE, NA, label)) %>% 
  mutate(nodelabs = na_if(nodelabs, "root"))  %>% 
  mutate(nodelabs = round(as.numeric(nodelabs), 2) * 100)

# add metadata
td <- treetib %>% left_join(sampdf %>% dplyr::rename(label = Species)) %>% 
  mutate(hyphae = "aaecd181-feb8-4203-8c64-f46384257e59",
         yeast = "e602729e-044f-4d2b-bab2-64a87a0b48c7")

# make figure
ffig.jacOTU <- ggtree(td, size = 1.5, ladderize = FALSE) +#, branch.length = "none") + 
  geom_rootedge(rootedge = 0.05, size = 1.5) +
  geom_tiplab(aes(image = hyphae, color = plotlab), geom = "phylopic",
              offset = -0.05, size = 0.09, alpha = 0.7)  +
  geom_tiplab(aes(image = yeast, color = plotlab), geom = "phylopic",
              offset = -0.05, size = 0.09) +
  geom_label(aes(label = nodelabs), nudge_x = -0.05, nudge_y = 0.2, fill = "white", label.size = 0, label.padding = unit(0, "lines")) + # WARNING IS OK
  geom_treescale(x = -0.1, y = 0, offset = 0.1) +
  
  
  # add stats
  geom_label(x = 0, y = 6, label = "Jaccard OTU \nnMC=0.53\nP=0.02", size = 5, label.size = NA) +
  #geom_tippoint(aes(color = plotlab, shape = Family), size = 7) +
  #theme_bw() + 
  scale_x_reverse(limits = c(1, -0.1)) + 
  theme(#legend.text = element_text(size = 16),#, color = "white"),
    #legend.title = element_text(size = 20),#, color = "white"),
    legend.position = "none",
    plot.margin = margin(l=0,r=0),
    plot.background = element_rect(fill = "white", color = NULL),
    panel.background = element_rect(fill = "white", color = NULL),
    legend.background = element_rect(fill = "white", color = NULL),
    legend.text = element_text(face = "italic", size = 12),
    legend.title = element_blank()) +
  scale_color_manual(values = ccol) 

## ---- fungal dendrogram: Bray Curtis, Genus ----

# read in
fd1 <- read.tree("data/topocong_input/bray_genusinR_rotated_meanrare1275.newick")
# fix tip labels
fd1$tip.label <- if_else(str_detect(fd1$tip.label, "troglody"), "pan_troglodytes", fd1$tip.label)
fd1$tip.label <- if_else(str_detect(fd1$tip.label, "procolobus"), "procolobus", fd1$tip.label)

# convert to treedata
treetib1 <- as_tibble(fd1) %>% as.treedata()

# wrangle
treetib1 <- treetib1 %>% 
  mutate(nodelabs = if_else(isTip == TRUE, NA, label)) %>% 
  mutate(nodelabs = na_if(nodelabs, "root"))  %>% 
  mutate(nodelabs = round(as.numeric(nodelabs), 2) * 100)

# add metadata
td1 <- treetib1 %>% left_join(sampdf %>% dplyr::rename(label = Species)) %>% 
  mutate(hyphae = "aaecd181-feb8-4203-8c64-f46384257e59",
         yeast = "e602729e-044f-4d2b-bab2-64a87a0b48c7")

# make figure
ffig.brayGenus <- ggtree(td1, size = 1.5, ladderize = FALSE) +#, branch.length = "none") + 
  geom_rootedge(rootedge = 0.05, size = 1.5) +
  geom_tiplab(aes(image = hyphae, color = plotlab), geom = "phylopic",
              offset = -0.05, size = 0.09, alpha = 0.7)  +
  geom_tiplab(aes(image = yeast, color = plotlab), geom = "phylopic",
              offset = -0.05, size = 0.09) +
  geom_label(aes(label = nodelabs), nudge_x = -0.05, nudge_y = 0.2, fill = "white", label.size = 0, label.padding = unit(0, "lines")) + # WARNING IS OK
  geom_treescale(x = -0.1, y = 0, offset = 0.1) +
  
  
  # add stats
  geom_label(x = 0, y = 5, label = "Bray-Curtis Genus \nnMC=0.71\nP=0.17", size = 5, label.size = NA) +
  #geom_tippoint(aes(color = plotlab, shape = Family), size = 7) +
  #theme_bw() + 
  scale_x_reverse(limits = c(1, -0.1)) + 
  theme(#legend.text = element_text(size = 16),#, color = "white"),
    #legend.title = element_text(size = 20),#, color = "white"),
    legend.position = "none",
    plot.margin = margin(l=0,r=0),
    plot.background = element_rect(fill = "white", color = NULL),
    panel.background = element_rect(fill = "white", color = NULL),
    legend.background = element_rect(fill = "white", color = NULL),
    legend.text = element_text(face = "italic", size = 12),
    legend.title = element_blank()) +
  scale_color_manual(values = ccol) 

## ---- fungal dendrogram: Jaccard, Genus ----

# read in
fd2 <- read.tree("data/topocong_input/jaccard_genusinR_rotated_meanrare1275.newick")
# fix tip labels
fd2$tip.label <- if_else(str_detect(fd2$tip.label, "troglody"), "pan_troglodytes", fd2$tip.label)
fd2$tip.label <- if_else(str_detect(fd2$tip.label, "procolobus"), "procolobus", fd2$tip.label)

# convert to treedata
treetib2 <- as_tibble(fd2) %>% as.treedata()

# wrangle
treetib2 <- treetib2 %>% 
  mutate(nodelabs = if_else(isTip == TRUE, NA, label)) %>% 
  mutate(nodelabs = na_if(nodelabs, "root"))  %>% 
  mutate(nodelabs = round(as.numeric(nodelabs), 2) * 100)

# add metadata
td2 <- treetib2 %>% left_join(sampdf %>% dplyr::rename(label = Species)) %>% 
  mutate(hyphae = "aaecd181-feb8-4203-8c64-f46384257e59",
         yeast = "e602729e-044f-4d2b-bab2-64a87a0b48c7")

# make figure
ffig.jacGenus <- ggtree(td2, size = 1.5, ladderize = FALSE) +#, branch.length = "none") + 
  geom_rootedge(rootedge = 0.05, size = 1.5) +
  geom_tiplab(aes(image = hyphae, color = plotlab), geom = "phylopic",
              offset = -0.05, size = 0.09, alpha = 0.7)  +
  geom_tiplab(aes(image = yeast, color = plotlab), geom = "phylopic",
              offset = -0.05, size = 0.09) +
  geom_label(aes(label = nodelabs), nudge_x = -0.05, nudge_y = 0.2, fill = "white", label.size = 0, label.padding = unit(0, "lines")) + # WARNING IS OK
  geom_treescale(x = -0.1, y = 0, offset = 0.1) +
  
  
  # add stats
  geom_label(x = 0, y = 6, label = "Jaccard Genus \nnMC=0.71\nP=0.17", size = 5, label.size = NA) +
  #geom_tippoint(aes(color = plotlab, shape = Family), size = 7) +
  #theme_bw() + 
  scale_x_reverse(limits = c(1, -0.1)) + 
  theme(#legend.text = element_text(size = 16),#, color = "white"),
    #legend.title = element_text(size = 20),#, color = "white"),
    legend.position = "none",
    plot.margin = margin(l=0,r=0),
    plot.background = element_rect(fill = "white", color = NULL),
    panel.background = element_rect(fill = "white", color = NULL),
    legend.background = element_rect(fill = "white", color = NULL),
    legend.text = element_text(face = "italic", size = 12),
    legend.title = element_blank()) +
  scale_color_manual(values = ccol) 

### ---- add together ----

### make one figure with the panels of each
leg <- get_legend(hostfig)

all <- ggarrange(hostfig.nolegend, ffig.jacOTU, hostfig.nolegend, ffig.brayGenus, hostfig.nolegend, ffig.jacGenus, ncol = 2, nrow = 3, widths = c(0.7, 1.1), labels = c("A.", "", "B.", "", "C.", "", "D.", ""), font.label = list(size = 14))

ggarrange(leg, all, ncol = 1, heights = c(0.1, 1)) +
  theme(plot.background = element_rect(fill = "white", color = NA))

ggsave(filename = "figures/supplementary_topocong.png", dpi = 300, height = 15, width = 12, units = "in")

