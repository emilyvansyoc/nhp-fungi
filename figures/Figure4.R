### sequence divergence figures 
#EVS 11/2024

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
library(phytools)
library(ggforce)
library(ggtree)
library(RColorBrewer)
library(deeptime)
library(ggnewscale)
library(TreeTools)
library(dplyr)
library(rphylopic)
library(treeio)


# get image file of sequence divergence analysis
load("data/image_sequencedivergence.RData")

# get colors
source("helpers/NHPColors.R")
newcol <- ccol[names(ccol) %in% c("H. sapiens", "P. troglodytes", "G. beringei", "G. gorilla")]
names(newcol) <- paste0("*", names(newcol), "*")
names(newcol) <- if_else(str_detect(names(newcol), "sapie"), "*H. sapien*", names(newcol))
names(newcol) <- if_else(str_detect(names(newcol), "troglo"), "*P. troglodyte*", names(newcol))


## load hominid tree
host <- read.tree("cophylogeny/cophylogeny_scripts/data/hominid_rotated.newick")

# remove mangabey
host1 <- drop.tip(host, tip = "Mangabey")

## ---- prettify fungal names and abbreviations ----


# order by Paco R2
load("cophylogeny/cophylogeny_scripts/data/final_sigs_otus_parafitpaco.RData")
## add letter abbreviations for plot
sigps <- sigps %>% 
  arrange(desc(pacor2)) %>% 
  mutate(letterab = LETTERS[1:nrow(sigps)+1])

#### get fungi that match timing
forplot <- clockest %>% 
  left_join(sigps %>% dplyr::select(otu, pacor2, letterab)) %>% 
  # order by decreasing r2
  arrange(desc(pacor2)) %>% 
  # make column of pretty-fied names that match codiv fig
  mutate(taxname = str_replace(Species, " Genus", " sp."),
         taxname = str_replace(taxname, " Order", " sp."),
         taxname = str_replace(taxname, "_", " ")) %>% 
  mutate(taxname = if_else( otu == "OTU_456", "Xylaria sp. 1", taxname),
         taxname = if_else(otu == "OTU_384", "Xylaria sp. 2", taxname)) %>% 
  
  mutate(timing = case_when(
    Label %in% "Between G. gorilla and H. sapiens" & between(mean, 7.1, 9.2) ~ "ontime",
    Label %in% "Between G. gorilla and H. sapiens" & mean < 7.1 ~ "early",
    Label %in% "Between G. gorilla and H. sapiens" & mean > 9.2 ~ "late"
  )) %>% 
  # shorten names of some fungi
  mutate(shortname = taxname) %>% 
  mutate(shortname = if_else(str_detect(shortname, "Saturnispora"), "S. diversa", shortname),
         shortname = if_else(str_detect(shortname, "hainan"), "N. hainanensis", shortname))

## ---- make scalar to convert hominid branchlengths to speciation time ----

# set the pan-homo split to 6
# set the panhomo-gorilla split to 8
# therefore, the scalar for BETWEEN the panhomo and gorilla is 2

# the branch length (x position) for homo/pan is 0.03593854
hp <- 0.03593854
scalar <- 2/hp

# create scale
sc <- data.frame(branchlength = host1$edge.length) %>% 
  mutate(toscale = branchlength*scalar) %>% 
  # set Gorilla split manually
  mutate(toscale = ifelse(between(branchlength, 0.06, 0.07), 7, toscale))



## change branch lengths of new tree
plottree <- host1
plottree$edge.length <- sc$toscale
#drop G beringei for clarity
plottree <- ape::drop.tip(plottree, tip = "Gorilla_beringei")

## ---- make dataframe for plotting mean times ----

## get dataframe for plotting the mean times
homgor <- forplot %>% 
  filter(Label == "Between G. gorilla and H. sapiens") %>% 
  mutate(y = 4.1, yend = 4.5) %>% 
  mutate(mean = if_else(mean > 10, 10.1, mean))

## get second dataframe for pan-gorilla
pangor <- forplot %>% 
  filter(Label == "Between G. gorilla and P. troglodytes") %>% 
  mutate(y = 1, yend = 1.4) %>% 
  mutate(mean = if_else(mean > 10, 10.1, mean))


## ----- function to make rectangles ----

# function to create a short dataframe for geom_rect
# courtesy of: https://www.sciencedirect.com/science/article/pii/S1055790322002330#s0135
interval_df_getter <- function(minage, maxage) {
  # bottom left, bottom right, top right, top left
  shape_data <- data.frame(x = c(-maxage, -minage, -minage, -maxage), y = c(-0.7, -0.7, 5.9, 5.9))
  return(shape_data)
}

## ---- get colors of fungal genera ----

### from: updated_figures_codiv.R
# use object: 'cols'

# get phylopic of a yeast
ypic <- "35db2572-a27d-4f83-9b45-11195d6fd5af"

# make dataframe
homgor <- homgor %>% 
  mutate(uuid = ypic) %>% 
  # make column for colors if they fall within the ranges
  mutate(rangecol = case_when(
    mean < 6 ~ "early",
    between(mean, 6, 7.09) ~ "rangeearly",
    between(mean, 7.1, 9.2) ~ "range",
    between(mean, 9.21, 10) ~ "rangelate",
    mean > 10 ~ "late"
  ))
pangor <- pangor %>% 
  mutate(uuid = ypic) %>% 
  # make column for colors if they fall within the ranges
  mutate(rangecol = case_when(
    mean < 6 ~ "early",
    between(mean, 6, 7.09) ~ "rangeearly",
    between(mean, 7.1, 9.2) ~ "range",
    between(mean, 9.21, 10) ~ "rangelate",
    mean > 10 ~ "late"
  ))


# make colors
#rangecols <- c("darkgray", "#FF781F", "#FF0000", "#FFBFB7", "darkgray")
rangecols <- c("darkgray", "lightblue", "darkblue", "lightblue", "darkgray")
names(rangecols) <- c("early", "rangeearly", "range", "rangelate", "late")

### 11/18; color by Figure 3
source("figures/Figure3_and_SuppFigureS3.R")
fcols <- data.frame(fcolor = unname(cols), forcolor = names(cols))
homgor <- homgor %>% 
  mutate(forcolor = if_else(str_detect(Genus, "Nigro"), "Xylaria/Nigrospora", Genus),
         forcolor = if_else(str_detect(forcolor, "Xylar"), "Xylaria/Nigrospora", forcolor),
         forcolor = if_else(str_detect(forcolor, "Pleo"), "Pleosporales", forcolor)) %>% 
  left_join(fcols)
pangor <- pangor %>% 
  mutate(forcolor = if_else(str_detect(Genus, "Nigro"), "Xylaria/Nigrospora", Genus),
         forcolor = if_else(str_detect(forcolor, "Xylar"), "Xylaria/Nigrospora", forcolor),
         forcolor = if_else(str_detect(forcolor, "Pleo"), "Pleosporales", forcolor)) %>% 
  left_join(fcols)

## ---- make plot ----


### SET UP
ggplot(plottree) +
  xlab(NULL) + ylab(NULL) +
  # set scales 
  scale_x_continuous(breaks = seq(-10, 0, 2),
                     minor_breaks = seq(-10, 0, 0.5),
                     labels = c(" > 10", abs(seq(-8, -2, 2)), "       0  Ma")) +
  # set y axis limits 
  ylim(c(-0.7, 6)) +
  ##### RECTANGLE SHAPES #####
# shape for "on time"
geom_shape(data = interval_df_getter(9.2, 7.1), aes(x = x, y = y), fill = alpha("#FF2400", 0.5), color = "black") +
  #annotate(geom = "text", x = -9.1, y= 5, size = 6, fontface = "bold", hjust = "left", label = "7.1 - 9.2 Ma") +
  # add text for homo-gorilla and pan-gorilla
  annotate(geom = "text", x = -10.1, y = 4.3, label = "Homo-Gorilla\nestimates", size = 6, fontface = "bold", hjust = "right") +
  annotate(geom = "text", x = -10.1, y = 1, label = "Pan-Gorilla\nestimates", size = 6, fontface = "bold", hjust = "right") +
  # shape for "close"
  geom_shape(data = interval_df_getter(9.23, 10), aes(x = x, y = y), fill = alpha("#FF781F", 0.3)) +
  geom_shape(data = interval_df_getter(6, 7.07), aes(x = x, y =y), fill = alpha("#FF781F", 0.3)) +
  # shape for "late"
  geom_shape(data = interval_df_getter(10.03, 11.8), aes(x = x, y = y), fill = alpha("grey", 0.3)) +
  # shape for "early"
  geom_shape(data = interval_df_getter(5.97, -0.5), aes(x = x, y = y), fill = alpha("grey", 0.3)) +
  #### ADD TREE ####
geom_tree(layout = "slanted", size = 2, position = position_nudge(x = -8, y = 0.8)) +
  # plot "root" by hand - geom_rootedge() does not have a 'position' argument
  geom_segment(x = -8, xend = -9, y = 2.55, linewidth = 2) +
  
  #### ADD STANDARD ERROR
  #geom_errorbar(data = homgor %>% filter(rangecol %in% c("rangeearly", "range", "rangelate")), aes(xmin = -mean-se, xmax = -mean+se, y = 4), width = 0.2) +
  ### ADD FUNGI PHYLOPICS ####
rphylopic::geom_phylopic(data = homgor %>% filter(rangecol %in% c("rangeearly", "range", "rangelate")), aes(uuid = ypic, x = -mean, y = 4, fill = forcolor), size = 0.3, color = "black") +
  rphylopic::geom_phylopic(data = pangor %>% filter(rangecol %in% c("rangeearly", "range", "rangelate")), aes(uuid = ypic, x = -mean, y = 1.5, fill = forcolor), size = 0.3, color = "black") +
  scale_fill_manual(values = cols) +
  ### ADD FUNGI NAMES BELOW ####
geom_text(data = homgor %>% filter(rangecol %in% c("rangeearly", "range", "rangelate")), aes(label = shortname, x = -mean, y = 4.2), angle = 270, size = 4.5, hjust = "right", fontface = "italic") +
  geom_text(data = pangor %>% filter(rangecol %in% c("rangeearly", "range", "rangelate")), aes(label = shortname, x = -mean, y = 1.3), angle = 270, size = 4.5, hjust = "left", fontface = "italic") +
  #### ADD HOMINID PHYLOPICS ####
# G beringei
#add_phylopic(uuid = "142e0571-3b5f-443d-a887-b572a224ea22", x = 0.2, y = 1, ysize = 0.6, fill = "#440154FF") +
# G gorilla
add_phylopic(uuid = "142e0571-3b5f-443d-a887-b572a224ea22", x = 0.2, y = 2, ysize = 0.6, fill = "#F781BF") +
  # human
  add_phylopic(uuid = "036b96de-4bca-408e-adb6-2154fcd724ef", x = -1.8, y = 3, ysize = 0.9, fill = "#984EA3") +
  # chimp
  add_phylopic(uuid = "7133ab33-cc79-4d7c-9656-48717359abb4", x = -2.8, y = 4, ysize = 0.6, fill = "#FF7F00") +
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(0.3, "cm"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.line.x = element_line(size = 1, linetype = "solid", color = "black"))

### save
ggsave(filename = "figures/Figure4.png", dpi = 300, height = 6, width = 12, units = "in")



