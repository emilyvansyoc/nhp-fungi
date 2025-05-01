### FIGURE 3

### SOME OF THIS IS PANELED IN ILLUSTRATOR and text colors of labels are changed in Illustrator (for now....)

library(tidyverse)
library(ape)
library(vegan)
library(ggtree)
library(ggpubr)
library(ggtext)
library(scales)
library(rphylopic)
library(treeio)

## load OTU trees
load("cophylogeny/cophylogeny_scripts/data/all_otu_subtrees.RData")

# get significant results
load("cophylogeny/cophylogeny_scripts/data/final_sigs_otus_parafitpaco.RData")

# get colors
source("helpers/NHPColors.R")
newcol <- ccol[names(ccol) %in% c("H. sapiens", "P. troglodytes", "G. beringei", "G. gorilla")]
names(newcol) <- paste0("*", names(newcol), "*")
names(newcol) <- if_else(str_detect(names(newcol), "sapie"), "*H. sapien*", names(newcol))
names(newcol) <- if_else(str_detect(names(newcol), "troglo"), "*P. troglodyte*", names(newcol))
# summary stats
sigps %>% arrange(desc(pacor2))

## ---- make function for making circular trees ---

myPlot <- function(mytree) {
  
  # break down hominid representation
  tipdf <- data.frame(tips = mytree$tip.label) %>% 
    mutate(tip1 = str_remove(tips, "_(\\d){1,4}$")) %>% 
    mutate(hominid = str_extract(tip1, "[:alpha:]{1,12}$")) %>% 
    column_to_rownames(var = "tips")
  tipdf %>% group_by(hominid) %>% count()
  
  # add sci name to tipdf
  tipdf <- tipdf %>% 
    mutate(italname = case_when(
      hominid == "Chimp" ~ "*P. troglodyte*",
      hominid == "Human" ~ "*H. sapien*",
      hominid == "LowGorilla" ~ "*G. gorilla*",
      hominid == "MountGorilla" ~ "*G. beringei*"
    )) 
  
  # simplify tip labels
  plottree <- mytree
  plottree$tip.label <- tipdf$italname
  td <- data.frame(rownames = plottree$tip.label,
                   hominid = plottree$tip.label)
  
  # plot circular
  tplot <- plottree %>% ggtree(ladderize = FALSE, layout = "fan", branch.length = "none") %<+% td + 
    geom_tippoint(aes(color = hominid), size = 3) +
    #scale_color_manual(values = italic.cols) +
    scale_color_manual(values = newcol) +
    theme(#legend.text = element_markdown(size = 18),
      #legend.title = element_blank(),
      legend.position = "none",
      plot.margin = margin(t=0, r=0, b=0, l=0, unit = "pt"),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA))
  
  # return plot
  return(tplot)
  
}

## ---- panel A: schematic ----

## load hominid tree
host <- read.tree("cophylogeny/cophylogeny_scripts/data/hominid_rotated.newick")

hostt <- drop.tip(host, tip = "Mangabey")
hostt$tip.label <- c("*G. gorilla*", "*G. beringei*", "*P. troglodyte*", "*H. sapien*" ) 
labdf <- data.frame(
  row.names = host$tip.label,
  name = host$tip.label,
  # uid = c("142e0571-3b5f-443d-a887-b572a224ea22", "142e0571-3b5f-443d-a887-b572a224ea22", ## GORILLA IS TWICE
  #     "036b96de-4bca-408e-adb6-2154fcd724ef", "7133ab33-cc79-4d7c-9656-48717359abb4"),
  forcolor = host$tip.label
) %>% 
  mutate(uid = case_when(
    name %in% "Gorilla_beringei" ~ "142e0571-3b5f-443d-a887-b572a224ea22",
    name %in% "Gorilla_gorilla" ~ "142e0571-3b5f-443d-a887-b572a224ea22",
    name %in% "Homo_sapien" ~ "036b96de-4bca-408e-adb6-2154fcd724ef",
    name %in% "P_troglodytes_schweinfurthii" ~ "7133ab33-cc79-4d7c-9656-48717359abb4",
    name %in% "Mangabey" ~ "eccbb404-c99f-41f9-8785-01a7f57f1269"
  ))

labdf1 <- labdf %>% filter(!name == "Mangabey") %>% mutate(label = case_when(
  name == "Gorilla_gorilla" ~ "*G. gorilla*",
  name == "Gorilla_beringei" ~ "*G. beringei*",
  name == "Homo_sapien" ~ "*H. sapien*",
  name == "P_troglodytes_schweinfurthii" ~ "*P. troglodyte*"
)) %>% mutate(dup = label) %>% rownames_to_column(var = "newname") %>% column_to_rownames(var = "dup") %>% dplyr::relocate(label)

hostt <- hostt %>% as_tibble() %>% as.treedata() %>% full_join(labdf1)

schem1 <- ggtree(hostt, size = 5, color = "black"#, branch.length = "none"
) + #%<+% labdf1 + ## the weird symbol is ggtree "attacher"
  geom_tiplab(aes(image = uid, color= label), geom = "phylopic", 
              offset = 0.05,
              size = c(0.09, # human 
                       0.14, 0.14, # gorillas
                       0.11
              )) +
  # show root
  geom_rootedge(rootedge = 0.05, size = 4) +
  xlim(-0.05,0.7) +
  scale_color_manual(values = newcol) +
  theme(legend.position = "none",
        plot.margin = margin(r=0))

#### make random tree to show codiversification schematic
# get a tree
set.seed(234)
#tr <- rtopology(n=12, rooted = TRUE)

#tr1 <- keep.tip(tr, tip = c("t4", "t3", "t24", "t11", "t10", "t8", "t21", "t15", "t2", "t9", "t1", "t18"))
#plot(tr1)
#tr1 <- tr
#tr$tip.label <- c("*G. beringei*", "*G. gorilla*", "*G. gorilla*", "*G. gorilla*", "*P. troglodyte*", "*H. sapien*", "*G. gorilla*", "*G. beringei*", "*G. beringei*", "*P. troglodyte*", "*H. sapien*", "*P. troglodyte*")
### save this random tree for re-use and tweaking the plot
load("cophylogeny/cophylogeny_scripts/data/image_for_codiv_schematic.RData")

# need unique tip labels
tr1$tip.label <- paste0("tip_", seq(1:12))


## get phylopic of a yeast
flabdf <- data.frame(
  row.names = tr1$tip.label,
  name = tr1$tip.label,
  uid = c("35db2572-a27d-4f83-9b45-11195d6fd5af"),
  forcolor = c("*G. gorilla*", "*G. beringei*", "*G. beringei*", "*G. beringei*", "*P. troglodyte*", "*H. sapien*", "*G. beringei*", "*G. gorilla*", "*G. gorilla*", "*P. troglodyte*", "*H. sapien*", "*P. troglodyte*")
)

### plot circle tree
### takes some tweaking to get the correct aspect ratios of our little yeast guys
df <- flabdf
df$svg <- lapply(df$uid, get_phylopic)

# make plot
schem2 <- ggtree(tr1, size = 5, layout = "circular", branch.length = "none", ladderize = FALSE) %<+% df +
  #geom_tiplab(aes(image = uid, color = forcolor)) +
  rphylopic::geom_phylopic(aes(img = svg, color = forcolor), size = 0.7) +
  scale_color_manual(values = newcol) +
  theme(#legend.text = element_markdown(size = 18),
    #legend.title = element_blank(),
    legend.position = "none",
    plot.margin = margin(t=0, r=0, b=0, l=0, unit = "pt"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA))

## ---- plot circular trees ----


# order by decreasing r2
sigps <- sigps %>% arrange(desc(pacor2))


# do by hand to get ordering right for now
### 11/13/24- update to keep only top 6 trees, move rest to supplementary
bigsig <- sigps[1:6,]
sigtrees <- alltrees[names(alltrees) %in% bigsig$otu]
sigtrees <- list(alltrees$OTU_169, alltrees$OTU_456, alltrees$OTU_329, alltrees$OTU_5023, alltrees$OTU_2339, alltrees$OTU_86)
names(sigtrees) <- bigsig$otu

# get "little sigs"
lsig <- sigps[7:11,]
lsigtrees <- alltrees[names(alltrees) %in% lsig$otu]
lsigtrees <- list(alltrees$OTU_559, alltrees$OTU_384, alltrees$OTU_282, alltrees$OTU_89, alltrees$OTU_12)

#orplot <- sigtrees[1:6] # plot all 10

plist <- lapply(sigtrees, myPlot)
names(plist) <- names(sigtrees)

### 10/31; add two random trees in white to make the spacing work
r1 <- rtree(5) %>% ggtree(ladderize = FALSE, layout = "fan", branch.length = "none", color = "white") + #%<+% td + 
  #geom_tippoint(aes(color = hominid), size = 2) +
  #scale_color_manual(values = newcol) +
  theme(#legend.text = element_markdown(size = 18),
    #legend.title = element_blank(),
    legend.position = "none",
    plot.margin = margin(t=0, r=0, b=0, l=0, unit = "pt"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA))
r2 <- rtree(5) %>% ggtree(ladderize = FALSE, layout = "fan", branch.length = "none", color = "white") + #%<+%# td + 
  #geom_tippoint(aes(color = hominid), size = 2) +
  #scale_color_manual(values = italic.cols) +
  theme(#legend.text = element_markdown(size = 18),
    #legend.title = element_blank(),
    legend.position = "none",
    plot.margin = margin(t=0, r=0, b=0, l=0, unit = "pt"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA))


# make a new plotlist
newplots <- c(list(r1), list(r2), plist)


ggarrange(plotlist = newplots, ncol = 4, nrow = 2)
#ggsave(filename = "updated_methods/fig_images/updated_codiv_otus_six_newcols.png", dpi = 600, height = 8, width = 12, units = "in")

## ---- Supplementary Figure 3 ----


# make tree points smaller for these trees (they are BIG)
myPlotSmall <- function(mytree) {
  
  # break down hominid representation
  tipdf <- data.frame(tips = mytree$tip.label) %>% 
    mutate(tip1 = str_remove(tips, "_(\\d){1,4}$")) %>% 
    mutate(hominid = str_extract(tip1, "[:alpha:]{1,12}$")) %>% 
    column_to_rownames(var = "tips")
  tipdf %>% group_by(hominid) %>% count()
  
  # add sci name to tipdf
  tipdf <- tipdf %>% 
    mutate(italname = case_when(
      hominid == "Chimp" ~ "*P. troglodyte*",
      hominid == "Human" ~ "*H. sapien*",
      hominid == "LowGorilla" ~ "*G. gorilla*",
      hominid == "MountGorilla" ~ "*G. beringei*"
    )) 
  
  # simplify tip labels
  plottree <- mytree
  plottree$tip.label <- tipdf$italname
  td <- data.frame(rownames = plottree$tip.label,
                   hominid = plottree$tip.label)
  
  # plot circular
  tplot <- plottree %>% ggtree(ladderize = FALSE, layout = "fan", branch.length = "none") %<+% td + 
    geom_tippoint(aes(color = hominid), size = 1) +
    #scale_color_manual(values = italic.cols) +
    scale_color_manual(values = newcol) +
    theme(#legend.text = element_markdown(size = 18),
      #legend.title = element_blank(),
      legend.position = "none",
      plot.margin = margin(t=0, r=0, b=0, l=0, unit = "pt"),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA))
  
  # return plot
  return(tplot)
  
}

#trees <- list(alltrees$OTU_79, alltrees$OTU_6)
#names(trees) <- c("OTU_79", "OTU_6")
plist1 <- lapply(lsigtrees, myPlotSmall)

# save together and make labels in Illustrator to color-code to main text figure
ggarrange(plotlist = plist1, ncol = 3, nrow = 2)
#ggsave(filename = "updated_methods/fig_images/codiv_supplementary_updated_newcol.png", dpi = 600, height = 6, width = 12, units = "in")

### ----- get full fungi tree and make schematic ----

ftree <- read.tree("figures/li2021_fungaltree/1672taxa_290genes_bb_1.treefile")


# get tip labels
ftips <- data.frame(tiplab = ftree$tip.label) %>% 
  mutate(genus = sapply(str_split(tiplab, "_"), `[`, 1))


# get the genera present in our dataset
load("cophylogeny/cophylogeny_scripts/data/all_genera_list.RData")

# subset
fsub <- ftips %>% 
  filter(genus %in% gen)

tr <- read.newick("figures/li2021_fungaltree/li2021_tree_mynodes.newick") %>% as.treedata()

td <- as_tibble(tr)

#### get Pleosporales nodes
# BLAST suggest these may belong to the Thyridariaceae family
tab <- readxl::read_xlsx(path = "figures/li2021_fungaltree/1-s2.0-S0960982221001391-mmc3.xlsx", sheet = "B")
# get Pleosporale
pleo <- tab %>% filter(order == "Pleosporales")
pleosub <- pleo %>% filter(tip_id %in% ftree$tip_id)
pleotree <- keep.tip(ftree, tip = pleo$tip_id)

# make dataframe of node labels by hand
# add labels to nodes
mynodes <- data.frame(
  genus = c("Aureobasidium", "Saturnispora", "Malassezia", 
            "Talaromyces", #"Geotrichum", 
            "Xylaria/Nigrospora", "Pleosporales", "Cladosporium",
            # 11/13 add
            "Lasiodiplodia"),
  mynode = c("node475", # aureobasidium
             "node247", # saturnispora
             "node414", # malassezia
             "node589", # talaro#"node327", 
             "node159", # xylaria/nigrospora
             # change pleosporales node to fit
             #"node449", 
             "node455", # pleo
             "node467", # cladosporium
             #"node448" falls in a weird spot on the tree - shift slightly
             "node463" # lasioplodia
  )) %>%  
  left_join(td, by = c("mynode" = "label"))



# get the default colors
library(scales)
cols <- hue_pal()(8)
names(cols) <- mynodes$genus


### build plot
allsigstree <- tr %>% ggtree() +
  geom_hilight(data = mynodes, aes(node = node, fill = genus),extend = 3.4, alpha = 0.9) +#,  type = "gradient", gradient.direction = "tr") +
  xlim(0, 3.4) +
  ggpubr::rotate() +
  theme(legend.position = "none") +
  scale_fill_manual(values = cols)
# save
#ggsave(filename = "updated_methods/fig_images/updated_ftree_elevensigs.png", dpi = 600)

##### only the fungi shown in figure 3
mynodes1 <- mynodes %>% filter(genus %in% c("Aureobasidium", "Xylaria/Nigrospora", "Talaromyces", "Pleosporales", "Saturnispora"))
somesigstree <- tr %>% ggtree() +
  geom_hilight(data = mynodes1, aes(node = node, fill = genus),extend = 3.4, alpha = 0.9) +#,  type = "gradient", gradient.direction = "tr") +
xlim(0, 3.4) +
  ggpubr::rotate() +
  theme(legend.position = "none") +
  scale_fill_manual(values = cols)
# save
#ggsave(filename = "updated_methods/fig_images/updated_ftree_elevensigs_forfig3.png", dpi = 600)
