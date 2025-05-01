### FIGURE 2

library(ggpubr)
library(ggtext)
library(rphylopic)
library(grid) # for rasterGrob
library(rphylopic)
library(ggh4x)
library(ggupset)
library(ggtext)
#library(patchwork)
library(ggplot2) #### see github issue; for now, downgrading ggplot2
# https://github.com/tidyverse/ggplot2/issues/5619

# source("R/colors.R")
# get updated colors
source("all_primate/R/NHPColors.R")

## get the analysis script
source("human_mycobiome/differential_relabundance.R")


## ---- panel A: upset plot ----

# get unique genera for each species
## get data
psgndf <- psg.named %>%
  tax_transform("pa") %>%
  ps_melt()

hum <- psgndf %>%
  group_by(Genus, ColName) %>%
  summarize(tot = sum(Abundance)) %>%
  ungroup() %>%
  mutate(pres = if_else(tot > 1, 1, 0)) %>%
  dplyr::select(-tot) %>%
  filter(pres == 1) %>%
  # add colors for text
  left_join(data.frame(ccol) %>% rownames_to_column(var = "ColName"), by = "ColName") %>%
  mutate(plotname = paste0("<i style='color:", ccol, "'>", ColName, "</b>"))

# show only certain intersections for clarity
myspecs <- sort(unique(hum$plotname))
myints <- list(
  myspecs[3],
  myspecs[7],
  myspecs[6],
  myspecs[4],
  myspecs[2],
  myspecs[8],
  myspecs[1],
  myspecs[5],
  c(myspecs[1], myspecs[2], myspecs[3], myspecs[4], myspecs[5], myspecs[6], myspecs[7], myspecs[8]),
  # just hominids
  c(myspecs[3], myspecs[4], myspecs[6], myspecs[7]),
  # just cercopithecoids
  c(myspecs[1], myspecs[2], myspecs[8])
)

# re-order the species
myspec.ord <- c("<i style='color:#984EA3'>H. sapiens</b>")

### BUILD PLOT ####
upplot <-  hum %>%
  group_by(Genus) %>%
  summarize(Species = list(plotname)) %>%
  ggplot(aes(x = Species)) +
  geom_bar(width = 0.6, fill = "darkgreen") +
  scale_x_upset( #intersections = myints, 
    #sets = myspecs, 
    n_intersections = 20
  ) +
  labs(x = NULL) +
  geom_text(stat = "count", aes(label = after_stat(count)), size = 6, vjust = -1) +
  scale_y_continuous(name = "Number of shared genera", limits = c(0, 60), breaks = c(0, 20, 40)) +
  theme_pubr() +
  axis_combmatrix(sep = "-") +
  theme_combmatrix(
    combmatrix.panel.striped_background.color.one = "grey85",
    combmatrix.panel.striped_background.color.two = "grey95",
    combmatrix.label.text = element_markdown(size = 16),
    combmatrix.label.make_space = FALSE,
    combmatrix.label.height = unit(200, "pt"),
    combmatrix.panel.point.size = 4
  ) +
  theme(text = element_text(size = 14),
        axis.title.y = element_text(size = 16, vjust = 2),
        plot.margin = margin(l=100, t= 40,  unit = "pt"))

### ----- panel B: dotplot ----

## make colors by hand
mycols <- c(ccol, "nosig" = "#D3D3D3")

#### Dotplot ####
p1 <- ggplot(lgm, aes(x = fc, y = ItalName, fill = issigcol, color = issigcol)) +
  facet_wrap(~gen.order,
             scales = "fixed", strip.position = "left", ncol = 1,
             dir = "v"
  ) +
  geom_segment(aes(x = 0, xend = fc, color = issigcol, yend = ItalName), position = position_dodge2(0.9), linetype = "solid", linewidth = 0.9) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0) +
  # rphylopic::geom_phylopic(aes(uuid = uuid), size = 2, position = position_dodge(0.9)) +
  scale_fill_manual(values = mycols, guide = "none") +
  scale_color_manual(values = mycols, guide = "none") +
  theme_pubr() +
  labs(x = "Log2 fold change") +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 16),
    panel.border = ggh4x::element_part_rect(side = "b", fill = NA, color = "black"),
    panel.spacing = unit(0.2, "cm"),
    strip.placement = "inside",
    strip.text.y.left = element_text(angle = 0, face = "italic", hjust = 1),
    text = element_text(size = 16),
    plot.margin = margin(t = 0, b = 0, unit = "pt")
  )

### make human phylopic and arrows on top of plot b 
h <- ggplot(lgm, aes(x = fc, y = 0, fill = issigcol)) +
  geom_blank() +
  geom_segment(aes(x = 0, xend = 5, y = -0.02, yend = -0.02), arrow = arrow(length = unit(0.2, "npc")), linewidth = 1) +
  geom_segment(aes(x = 0, xend = -5, y = -0.02, yend = -0.02), arrow = arrow(length = unit(0.2, "npc")), linewidth = 1) +
  geom_text(x = 3, y = 0.01, label = "More abundant in \nH. sapiens", fontface = "italic", size = 5, check_overlap = TRUE, lineheight = 0.8) +
  geom_text(x = -3, y = 0.01, label = "Less abundant in \nH. sapiens", fontface = "italic", size = 5, check_overlap = TRUE, lineheight = 0.8) +
  add_phylopic(uuid = "036b96de-4bca-408e-adb6-2154fcd724ef", x = 0, y = 0, fill = "#984EA3", ysize = 0.1) +
  theme_void() +
  ylim(c(-0.025, 0.025)) +
  theme(plot.margin = margin(t = 0.5, b = 0, unit = "cm"))


### ---- panel C: tiles ----

# get top prev by total prev
topp <- psg.named %>% tax_top(n = 20, by = "prev", rank = "Genus")
patop <- psgndf %>% filter(OTU %in% topp)


# Create heatmap of OTU abundance by Species
# First get the order of OTUs based on prevalence in humans
otu_order <- patop %>%
  # filter(Species == "homo_sapiens") %>%
  group_by(OTU) %>%
  summarize(all_prev = sum(Abundance)) %>% # or mean(Abundance) depending on your data
  arrange(desc(all_prev)) %>%
  slice_head(n = 16) %>%
  pull(OTU)

# set colors
fillcols <- c("#ECECEC", "#191970")
names(fillcols) <- c(0, 1)

# Create the heatmap
matplot <- patop %>%
  mutate(Species = case_when(
    Species %in% "homo_sapiens" ~ "H. sapiens",
    Species %in% "pan_troglodytes" ~ "P. troglodytes",
    Species %in% "gorilla_gorilla" ~ "G. gorilla",
    Species %in% "gorilla_beringei" ~ "G. beringei",
    Species %in% "cercocebus_agilis" ~ "C. agilis",
    Species %in% "indri_indri" ~ "I. indri",
    Species %in% "papio_cynocephalus" ~ "P. cynocephalus",
    Species %in% "procolobus" ~ "P. gordonum"
  )) %>%
  filter(OTU %in% otu_order) %>%
  group_by(OTU, Species) %>%
  arrange(desc(Abundance)) %>%
  mutate(
    OTU = factor(OTU, levels = rev(otu_order)),
    Species = factor(Species, ordered = TRUE, levels = c("H. sapiens", "P. troglodytes", "G. gorilla", "G. beringei", "P. gordonum", "C. agilis",  "P. cynocephalus", "I. indri" ))
  ) %>% # reorder OTUs
  ggplot(aes(x = Sample, y = OTU, fill = factor(Abundance))) +
  geom_tile() +
  ggh4x::facet_wrap2(~Species,
                     scales = "free_x", ncol = 8, nrow = 1,
                     strip = strip_themed(background_x = elem_list_rect(fill = c("#984EA3", "#FF7F00", "#F781BF", "#A65628", "#4DAF4A", "#FFDE21", "#377EB8", "#E41A1C")))
  ) +
  scale_fill_manual(values = fillcols) +
  # theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(face = "italic", size = 14),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.01, "cm"),
    strip.text = element_text(color = "white", face = "italic", size = 14),
    axis.title = element_text(size = 14),
    panel.grid = element_blank(),
    legend.position = "none" # removes legend
  ) +
  labs(
    x = "Individual",
    y = "Prevalent genera"
  )

### ---- combine and save ----

### there are some dependency issues with patchwork; call namespace for just these two functions to avoid loading
patch1 <- patchwork::free(h, side = "tb", type = c("space")) / p1 +
  patchwork::plot_layout(heights = c(0.1, 1))

# arrange

top <- ggarrange(upplot, patch1, ncol = 2, labels = c("A.", "B."), font.label = list(size = 16))
ggarrange(top, matplot, ncol = 1, heights = c(1, 0.8), labels = c("", "C."), font.label = list(size = 16)) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave(filename = "figures/figure2.png", dpi = 300, height = 13, width = 16, units = "in")

