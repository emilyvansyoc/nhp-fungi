## supplementary figure 2

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

# calculate log2 change
lg1 <- psg.named %>%
  tax_select(unique(sig_mods$taxon), strict_matches = TRUE) %>% 
  tax_select(sig_taxa, strict_matches = TRUE, deselect = TRUE) %>%
  tax_transform("compositional") %>%
  tax_transform("log2", zero_replace = "halfmin", chain = TRUE) %>%
  ps_melt() # %>%
# left_join(all1 %>% dplyr::select(taxon, comp, padj, estimate), by = c("OTU" = "taxon", "SciName" = "comp")) %>%
# mutate(is.sig = if_else(padj < 0.05, TRUE, FALSE))

# get the means
lgm1 <- lg1 %>%
  group_by(OTU, Species) %>%
  summarize(avg = mean(Abundance)) %>%
  ungroup() %>%
  pivot_wider(names_from = Species, values_from = avg) %>%
  # to calculate log2 fold change, want log2(A) - log2(B)
  mutate(
    GB = homo_sapiens - gorilla_beringei,
    GG = homo_sapiens - gorilla_gorilla,
    PT = homo_sapiens - pan_troglodytes,
    CA = homo_sapiens - cercocebus_agilis,
    II = homo_sapiens - indri_indri,
    PC = homo_sapiens - papio_cynocephalus,
    PR = homo_sapiens - procolobus
  ) %>%
  dplyr::select(OTU, G_beringei = GB, G_gorilla = GG, P_troglodytes = PT, C_agilis = CA, I_indri = II, P_cynocephalus = PC, P_colobus = PR) %>%
  pivot_longer(!OTU, names_to = "SciName", values_to = "fc") %>%
  mutate(comp = paste(OTU, SciName)) %>%
  left_join(all_mods, by = "comp") %>%
  drop_na() %>% 
  mutate(issig = if_else(padj < 0.05, TRUE, FALSE)) %>%
  mutate(ItalName = case_when(
    SciName %in% "G_gorilla" ~ "G. gorilla",
    SciName %in% "G_beringei" ~ "G. beringei",
    SciName %in% "P_troglodytes" ~ "P. troglodytes",
    SciName %in% "C_agilis" ~ "C. agilis",
    SciName %in% "I_indri" ~ "I. indri",
    SciName %in% "P_cynocephalus" ~ "P. cynocephalus",
    SciName %in% "P_colobus" ~ "P. colobus"
  )) %>%
  mutate(SciName = factor(SciName, ordered = TRUE, levels = c("P_troglodytes", "G_gorilla", "G_beringei", "C_agilis", "I_indri", "P_cynocephalus", "P_colobus"))) %>%
  # add gray color for insignificant comparisons
  mutate(issigcol = if_else(issig == TRUE, ItalName, "nosig"))

lgm1 <- lgm1 %>%
  # order fungal genera
  mutate(gen.order = factor(OTU, ordered = TRUE, levels = lgm1 %>% arrange(desc(fc)) %>% distinct(OTU) %>% pull(OTU)))

### plot


#### Dotplot ####
supp <- ggplot(lgm1, aes(x = fc, y = ItalName, fill = issigcol, color = issigcol)) +
  facet_wrap(~gen.order,
             scales = "free_y", strip.position = "left", ncol = 3,
             dir = "v"
  ) +
  geom_segment(aes(x = 0, xend = fc, color = issigcol, yend = ItalName), position = position_dodge2(0.9), linetype = "solid", linewidth = 0.9) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0) +
  # rphylopic::geom_phylopic(aes(uuid = uuid), size = 2, position = position_dodge(0.9)) +
  scale_fill_manual(values = mycols, guide = "none") +
  scale_color_manual(values = mycols, breaks = c("P. troglodytes", "P. cynocephalus",  "I. indri","P. gordonum",  "G. gorilla", "G. beringei", "C. agilis"), name = "NHP") +
  theme_pubr() +
  labs(x = "Log2 fold change") +
  theme(
    legend.position = "right",
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    panel.border = ggh4x::element_part_rect(side = "b", fill = NA, color = "black"),
    panel.spacing = unit(0.2, "cm"),
    strip.placement = "inside",
    strip.text.y.left = element_text(angle = 0, face = "italic", hjust = 1),
    text = element_text(size = 12),
    plot.margin = margin(t = 0, b = 0, unit = "pt")
  )

## add human
h1 <- ggplot(lgm1, aes(x = fc, y = 0, fill = issigcol)) +
  geom_blank() +
  geom_segment(aes(x = 0, xend = 5, y = -0.02, yend = -0.02), arrow = arrow(length = unit(0.1, "npc")), linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = -5, y = -0.02, yend = -0.02), arrow = arrow(length = unit(0.1, "npc")), linewidth = 0.5) +
  geom_text(x = 3, y = 0.01, label = "More abundant in \nH. sapiens", fontface = "italic", size = 3, check_overlap = TRUE, lineheight = 0.8) +
  geom_text(x = -3, y = 0.01, label = "Less abundant in \nH. sapiens", fontface = "italic", size = 3, check_overlap = TRUE, lineheight = 0.8) +
  add_phylopic(uuid = "036b96de-4bca-408e-adb6-2154fcd724ef", x = 0, y = 0, fill = "#984EA3", ysize = 0.1) +
  theme_void() +
  ylim(c(-0.025, 0.025)) +
  theme(plot.margin = margin(t = 1, b = 0, unit = "cm"),
        plot.background = element_rect(fill = "white", color = "white"))

#patch1 <- free(h1, side = "b", type = c("space")) / supp +
# plot_layout(heights = c(0.05, 1), axes = "collect") +
#plot_annotation(theme = theme(plot.margin = margin(t = 0.5, unit = "cm")))

par1 <- ggarrange(h1, supp, ncol = 1, heights = c(0.08, 1))
par1 <- par1 + theme(plot.background = element_rect(fill = "white", color = "white"),
                     plot.margin = margin(t = 1, b = 0, unit = "cm"))

ggsave(par1, filename = "figures/SupplementaryFigure2.png", dpi = 300, height = 17, width = 17, units = "in")
