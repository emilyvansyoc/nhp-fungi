#### generating Figure 1 plot ####
# EVS 4/2025


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
library(ape)
library(ggtree)
library(treeio)
library(gridExtra)
library(rphylopic)
library(ggimage)
library(ggpubr)

# set colors
library(RColorBrewer)
library(scales)


# get colors
source("helpers/NHPColors.R")

### source several previous scripts to get RData objects to build plots
source("pairwise_braycurtis/pairwise_distances.R")

### ---- PANEL A: RELATIVE ABUNDANCE BARPLOT ----

# pretty-fy names
filt99 <- filt99 %>% 
  ps_mutate(SciName = case_when(
    Species == "cercocebus_agilis" ~ "C. agilis",
    Species == "pan_troglodytes" ~ "P. troglodytes",
    Species == "gorilla_gorilla" ~ "G. gorilla",
    Species == "homo_sapiens" ~ "H. sapiens",
    Species == "gorilla_beringei" ~ "G. beringei",
    Species == "indri_indri" ~ "I. indri",
    Species == "papio_cynocephalus" ~ "P. cynocephalus",
    Species == "procolobus" ~ "P. gordonum"
  ))

## plot with hierarchical colors
# this code is mESSY BUT NOT MINE: https://david-barnett.github.io/microViz/articles/web-only/compositions.html


hueRank <- "Class"
hueRankPlural <- "Class"
shadeRank <- "Genus"

# Sort phyloseq at lower, and then higher ranks
pseq2 <- filt99 %>%
  tax_fix() %>% 
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" families will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(pseq2)[, hueRank]),
  shade = as.vector(tt_get(pseq2)[, shadeRank]),
  counts = taxa_sums(otu_get(pseq2))
)

hierarchicalPalInfo <- hierarchicalPalInfo %>%
  dplyr::mutate(
    hue = forcats::fct_other(
      f = hue, keep = unique(hue)[seq_len(nHues)],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
  ) %>%
  dplyr::group_by(hue) %>%
  dplyr::mutate(
    shade = forcats::fct_other(
      f = shade, keep = unique(shade)[seq_len(nShades - 1)],
      other_level = "Other"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
    Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
  )
## remove "others" by hand
hierarchicalPalInfo <- hierarchicalPalInfo %>% 
  mutate(hue = if_else(hue == "Other Class", "Other", hue),
         Taxa = if_else(str_detect(Taxa, "Other Class"), "Other genera", Taxa)) 

hierarchicalPalMatrix <- matrix(
  data = sapply(
    X = seq(from = 30, to = 75, length.out = nShades),
    FUN = function(l) scales::hue_pal(l = l, h.start = 30)(n = nHues)
  ),
  byrow = TRUE, ncol = nHues
)
hierarchicalPalMatrix <- cbind(hierarchicalPalMatrix, grey.colors(n = nShades))

hierarchicalPal <- hierarchicalPalMatrix %>%
  as.vector() %>%
  setNames(unique(hierarchicalPalInfo$Taxa)) 
#hierarchicalPal <- hierarchicalPal[!is.na(hierarchicalPal)]
#hierarchicalPal <- hierarchicalPal[!is.na(names(hierarchicalPal))]




hierarchicalPal <- hierarchicalPal[!is.na(names(hierarchicalPal))]
hierarchicalPal[names(hierarchicalPal) == "Other genera"] <- "#EEEEEE"
hierarchicalPal <- hierarchicalPal[!is.na(names(hierarchicalPal))]
names(hierarchicalPal) <- str_trim(names(hierarchicalPal))

tax_palette_plot(hierarchicalPal) +
  theme(axis.text.y.left = element_text(family = "Arial"))

bplot <- pseq2 %>%
  # order individuals by composition (roughly)
  ps_seriate(rank = "Genus") %>%
  ps_mutate(subject = sample_names(pseq2),
            subject = factor(subject, levels = unique(subject))) %>% 
  # order species by phylogeny
  ps_mutate(SciName = factor(SciName, ordered = TRUE, levels = c("G. beringei", "G. gorilla", "P. troglodytes", "H. sapiens", "P. gordonum", "C. agilis", "P. cynocephalus", "I. indri"))) %>% 
  ps_get() %>%
  tax_mutate("Class: Genus" = str_trim(hierarchicalPalInfo$Taxa), .keep = "none") %>%
  comp_barplot(
    tax_level = "Class: Genus", n_taxa = length(hierarchicalPal),
    tax_order = "asis", #palette = hierarchicalPal, 
    bar_width = 0.975, bar_outline_width = 0.01,
    label = NULL, x = "subject"
  ) +
  facet_wrap2(~SciName, ncol = 8, nrow = 1, scales = "free_x", 
              strip = strip_themed(background_x = elem_list_rect(fill = c("#A65628", "#F781BF", "#FF7F00", "#984EA3", "#4DAF4A", "#FFDE21", "#377EB8", "#E41A1C")))) +
  labs(x = "Individual", y = "Relative Abundance") +
  scale_fill_manual(values = hierarchicalPal, breaks = names(hierarchicalPal)[1:12], guide = guide_legend(ncol = 3)) +
  #guides(fill = guide_legend(ncol = 4, label.theme = element_text(size = 11), label.hjust = 0)) +
  #coord_flip() +
  ggplot2::theme(legend.text = element_text(size = 12),
                 text = element_text(size = 14),
                 legend.position = "bottom",
                 #legend.location = "plot",
                 legend.title = element_blank(),
                 legend.box.margin = margin(t=0,b=0, unit = "cm"),
                 legend.box.spacing = unit(0.01, units = "cm"),
                 #legend.key.spacing = unit(0.05, units = "cm"),
                 panel.spacing = unit(0.01, units = "cm"),
                 #strip.background = element_rect(fill = "white", color = "black"),
                 strip.text = element_text(face = "bold.italic", color = "white"),
                 axis.ticks.x = element_blank(),
                 axis.title.x = element_blank(),
                 legend.key.size = unit(0.7, "lines")) 

### get legend
bleg <- get_legend(bplot)

### plot without legend
bplot.nolegend <- bplot + theme(legend.position = "none")

#### ---- PANEL B: TOPO CONG ----

## setup
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


#### get phylopics
h.img <- get_phylopic(uuid = "036b96de-4bca-408e-adb6-2154fcd724ef", preview = TRUE)
g.img <- get_phylopic(uuid = "142e0571-3b5f-443d-a887-b572a224ea22", preview = TRUE)
c.img <- get_phylopic(uuid = "7133ab33-cc79-4d7c-9656-48717359abb4", preview = TRUE)
f.img <- get_phylopic(uuid = "aaecd181-feb8-4203-8c64-f46384257e59", preview = TRUE) # hyphae
f.img1 <- get_phylopic(uuid = "e602729e-044f-4d2b-bab2-64a87a0b48c7", preview = TRUE) # yeast
pap.img <- get_phylopic(uuid = "72f2f854-f3cd-4666-887c-35d5c256ab0f", preview = TRUE)
proc.img <- get_phylopic(uuid = "525fde4f-596e-414f-8a2c-0837aa75792c", preview = TRUE)
ind.img <- get_phylopic(uuid = "b28bef2c-0f7a-45e7-a815-f67cd1a98d8f", preview = TRUE)
m.img <- get_phylopic(uuid = "eccbb404-c99f-41f9-8785-01a7f57f1269", preview = TRUE)

#### HOST TREE
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
  # add stats
  geom_label(x = 10, y = 7, label = "nMC=0.53\nP=0.02", size = 5, label.size = NA) +
  geom_treescale(x = 0, y = 0, offset = 0.1) + 
  xlim(-10, 80) +
  theme(#legend.text = element_text(size = 16),
    #legend.title = element_text(size = 20),
    legend.position = "bottom",
    plot.margin = margin(l = 0,r=0),
    plot.background = element_rect(fill = "white", color = NULL),
    panel.background = element_rect(fill = "white", color = NULL),
    legend.background = element_rect(fill = "white", color = NULL),
    legend.text = element_text(face = "italic", size = 12),
    legend.title = element_blank()) +
  scale_color_manual(values = ccol) +
  scale_size_manual(values = csize)

#### FUNGAL DENDROGRAM 
# read in
fd <- read.tree("data/topocong_input/bray_OTU_rotated_meanrare1837.newick")
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
td <- treetib %>% left_join(sampdf %>% dplyr::rename(label = Species))
# add fungi and yeast phylopics
td <- td %>% 
  mutate(hyphae = "aaecd181-feb8-4203-8c64-f46384257e59",
         yeast = "e602729e-044f-4d2b-bab2-64a87a0b48c7")

# make figure
ffig <- ggtree(td, size = 1.5, ladderize = FALSE) +#, branch.length = "none") + 
  geom_rootedge(rootedge = 0.05, size = 1.5) +
  geom_tiplab(aes(image = hyphae, color = plotlab), geom = "phylopic",
              offset = -0.05, size = 0.09, alpha = 0.7)  +
  geom_tiplab(aes(image = yeast, color = plotlab), geom = "phylopic",
              offset = -0.05, size = 0.09) +
  geom_label(aes(label = nodelabs), nudge_x = -0.05, nudge_y = 0.2, fill = "white", label.size = 0, label.padding = unit(0, "lines")) + # WARNING IS OK
  geom_treescale(x = -0.1, y = 0, offset = 0.1) +
  #geom_tippoint(aes(color = plotlab, shape = Family), size = 7) +
  #theme_bw() + 
  scale_x_reverse(limits = c(1, -0.1)) + 
  theme(#legend.text = element_text(size = 16),#, color = "white"),
    #legend.title = element_text(size = 20),#, color = "white"),
    legend.position = "bottom",
    plot.margin = margin(l=0,r=0),
    plot.background = element_rect(fill = "white", color = NULL),
    panel.background = element_rect(fill = "white", color = NULL),
    legend.background = element_rect(fill = "white", color = NULL),
    legend.text = element_text(face = "italic", size = 12),
    legend.title = element_blank()) +
  scale_color_manual(values = ccol) 

### add both
pboth <- ggarrange(hostfig, ffig, ncol = 2, common.legend = TRUE, legend = "none", widths = c(1, 1))

#### ----- PANEL C: INTRASPECIFIC ----


# make color vector
sci.cols <- ccol
names(sci.cols) <- paste0("Within-", names(sci.cols))

# add "between" colors (all white)
f.cols <- c(sci.cols, rep("white", length(sci.cols)))
names(f.cols)[9:16] <- paste0("Between-", names(ccol))


# get stats for plot #### since all are signiicant above, use t test for plotting simplicity
stats <- dat %>%
  group_by(plotcol) %>%
  # anova_test(formula = Distance ~ is.within) %>%
  rstatix::t_test(Distance ~ is.within) %>%
  # wonky code to make it only have one asterisk for simplicity
  add_significance(cutpoints = c(0, 0.05, 1, 2), symbols = c("*", "+", "ns")) %>%
  add_xy_position() %>%
  mutate(y.position = 0.15) %>% 
  mutate(plotcol = factor(plotcol, ordered = TRUE, levels = levs))

##### MAKE PLOT
plota <- ggplot(data = dat, aes(x = is.within, y = Distance)) +
  geom_jitter(
    width = 0.2, color = "darkgrey", alpha = 0.3
  ) +
  geom_boxplot(#aes(fill = is.within),
    aes(fill = fulllab),
    size = 1.5, outlier.shape = NA, notch = FALSE, alpha = 0.7
  ) +
  coord_flip() +
  ggh4x::facet_wrap2(~plotcol, # strip = strip_themed(background_y = elem_list_rect(fill = c("#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF"))),
                     ncol = 1, strip.position = "left"
  ) +
  # add phylopics
  rphylopic::geom_phylopic(data = ppic, aes(uuid = uids, color = plotcol), x = 1.5, y = 0.1, size = 2) +
  ylim(c(0.05, 1.01)) +
  labs(y = "Bray-Curtis pairwise distances", x = "") +
  stat_pvalue_manual(data = stats %>% mutate(y.position = 0.2), label = "p.signif", hide.ns = TRUE, tip.length = 0, size = 10, coord.flip = TRUE) + # brackets
  scale_fill_manual(values = f.cols) +
  # ylim(c(0, 1.15)) +
  scale_x_discrete(position = "top") +
  theme_pubr(base_size = 14) +
  theme( # strip.background = element_rect(fill = c("#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF", alpha = 0.9), color = NA),
    #strip.text = element_text(size = 16, face = "italic", color = "black"),
    strip.text = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.position = "none",
    panel.spacing = unit(0.05, "cm")
  ) +
  scale_color_manual(values = ccol)

## ---- PANEL D: INTERSPECIFIC ----

# change label for test; easier to see
mod <- aov(Distance ~ labplot, data = humv)
t <- TukeyHSD(mod)
# use multcompview to get letters easily
cld <- multcompLetters4(mod, t, reversed = TRUE)
cldf <- as.data.frame.list(cld$labplot) %>%
  rownames_to_column(var = "labplot") %>%
  dplyr::select(labplot, Letters) %>%
  # get positions
  full_join(humv %>% group_by(labplot) %>% get_summary_stats(Distance, type = "common") %>% dplyr::select(median, labplot)) %>% 
  arrange(desc(median))


## make colors letters
coldf <- ccol %>% as.data.frame() %>% rownames_to_column(var = "labplot")
names(coldf) <- c("labplot", "col")
humv1 <- humv %>% 
  full_join(coldf %>% filter(!labplot == "H. sapiens")) %>% 
  mutate(labplot = factor(labplot, ordered = TRUE, levels = cldf$labplot)) %>% 
  #mutate(fulllab = paste0("Between <b style='color:#984EA3'>*H. sapiens*</b> and <b style='color:", col, "'>*", labplot, "*</b>"))
  mutate(fulllab = paste0("<b style='color:#984EA3'>*H. sapiens*</b> vs. <b style='color:", col, "'>*", labplot, "*</b>")) 
cldf1 <- cldf %>% left_join(humv1 %>% dplyr::select(fulllab, labplot) %>% distinct()) %>% arrange(desc(Letters))
humv1 <- humv1 %>% mutate(fulllab = factor(fulllab, ordered = TRUE, levels = cldf1$fulllab))

#### plot
plotb <- ggplot(data = humv1, aes(x = fct_reorder(fulllab, desc(Distance)), y = Distance)) +
  
  geom_point(position = position_jitter(width = 0.2, seed = 123), alpha = 0.5, color = "darkgrey") +
  geom_boxplot(size = 1.5, outlier.shape = NA, notch = FALSE, fill = "lightcyan4", alpha = 0.7) +
  geom_text(data = cldf1, aes(x = fulllab, label = Letters, y = 1.05),
            size = 6) +
  #ylim(-0.1, 1.1) +
  scale_y_continuous(limits = c(0.5, 1.1)) +
  coord_flip() +
  theme_pubr(base_size = 12) +
  labs(x = "", y = "Bray-Curtis pairwise distances") +
  theme(axis.text.y = element_markdown())

### change axis of barplot to make it easier to see the small spread
plotb1 <- ggplot(data = humv1 %>% 
                   mutate(newdis = if_else(Distance <= 0.8, 0.8, Distance)), aes(x = fct_reorder(fulllab, desc(Distance)), y = newdis)) +
  
  geom_point(position = position_jitter(width = 0.2, seed = 123), alpha = 0.5, color = "darkgrey") +
  geom_boxplot(size = 1.5, outlier.shape = NA, notch = FALSE, fill = "lightcyan4", alpha = 0.7) +
  geom_text(data = cldf1, aes(x = fulllab, label = Letters, y = 1.02),
            size = 6) +
  #ylim(-0.1, 1.1) +
  scale_y_continuous(limits = c(0.8, 1.05), breaks = c(0.8, 0.9, 1), labels = c("<0.8", "0.9", "1")) +
  coord_flip() +
  theme_pubr(base_size = 12) +
  labs(x = "", y = "Bray-Curtis pairwise distances") +
  theme(axis.text.y = element_markdown())

### ---- ADD TOGETHER AND SAVE ----

top1 <- ggarrange(bplot, pboth, ncol = 2, labels = c("A.", "B."), widths = c(1.5, 1), font.label = list(size = 14))
bottom <- ggarrange(plota, plotb1, labels = c("C.", "D."), ncol = 2, #widths = c(0.7, 1), 
                    font.label = list(size = 14))
all1 <- ggarrange(top1, bottom, ncol = 1, nrow = 2, heights = c(1, 0.8))
ggsave(all1, filename = "figures/figure1.png", dpi = 300, height = 11, width = 18, units = "in")
