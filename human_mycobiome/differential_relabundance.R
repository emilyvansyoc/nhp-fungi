## differential relative abundance: human-v-NHP
# EVS 4/2025

library(microViz)
library(phyloseq)
library(tidyverse)


## get data
load("data/phylo_OTU99_ITS2.RData")

# rename to scientific name so tip labels match
psfilt <- filt99 %>%
  ps_mutate(ishuman = if_else(str_detect(Species, "sapien"), 1, 0)) %>%
  ps_mutate(ColName = case_when(
    Species %in% "homo_sapiens" ~ "H. sapiens",
    Species %in% "pan_troglodytes" ~ "P. troglodytes",
    Species %in% "gorilla_gorilla" ~ "G. gorilla",
    Species %in% "gorilla_beringei" ~ "G. beringei",
    Species %in% "cercocebus_agilis" ~ "C. agilis",
    Species %in% "indri_indri" ~ "I. indri",
    Species %in% "papio_cynocephalus" ~ "P. cynocephalus",
    Species %in% "procolobus" ~ "P. gordonum"
  ))

### get only named genera for rel abund testing
gen <- get_taxa_unique(psfilt, "Genus")
gen <- gen[!is.na(gen)]
psg.named <- psfilt %>%
  tax_select(gen, ranks_searched = "Genus", strict_matches = TRUE) %>%
  tax_agg("Genus") # 515 genera


## --- what are the core prevalent genera? ----

## what is the prevalence of the top taxa across species?
# calculate prevalence within each species
df <- psg.named %>%
  tax_transform("pa") %>%
  ps_melt() %>%
  dplyr::select(OTU, Sample, Abundance, Species, ishuman)
prev <- df %>%
  group_by(OTU, Species) %>%
  summarize(tot = sum(Abundance))
# get species totals to calculate prevalence
stot <- psfilt %>%
  samdat_tbl() %>%
  group_by(Species) %>%
  count()
# calculate prevalence
prevg <- prev %>%
  mutate(prevalence = case_when(
    Species == "gorilla_beringei" ~ tot / stot$n[stot$Species == "gorilla_beringei"],
    Species == "gorilla_gorilla" ~ tot / stot$n[stot$Species == "gorilla_gorilla"],
    Species == "pan_troglodytes" ~ tot / stot$n[stot$Species == "pan_troglodytes"],
    Species == "homo_sapiens" ~ tot / stot$n[stot$Species == "homo_sapiens"],
    Species == "cercocebus_agilis" ~ tot / stot$n[stot$Species == "cercocebus_agilis"],
    Species == "indri_indri" ~ tot / stot$n[stot$Species == "indri_indri"],
    Species == "papio_cynocephalus" ~ tot / stot$n[stot$Species == "papio_cynocephalus"],
    Species == "procolobus" ~ tot / stot$n[stot$Species == "procolobus"],
  )) %>%
  dplyr::select(-tot) %>%
  filter(prevalence > 0) # %>%
# group_by(OTU) %>%
# count()  %>% arrange(desc(n)) #%>%
# filter(n == 8) # only two are present in all 8 species at 10% prevalence; 13 are present in all 8 species at any prevalence

## ---- get rel abund for human-specific taxa ----

psgndf <- psg.named %>%
  tax_transform("pa") %>%
  ps_melt()

# how many genera shared across all host species?
shared.g <- psgndf %>%
  group_by(Genus, Species) %>%
  summarize(tot = sum(Abundance)) %>%
  ungroup() %>%
  mutate(pres = if_else(tot > 1, 1, 0)) %>%
  group_by(Genus) %>%
  summarize(prev = sum(pres)) %>%
  arrange(desc(prev)) # 4 shared; Cladosporium, Ganoderma, Penicillium, Talaromyces
# how many are unique to a species?
length(which(shared.g$prev == 1)) # 221

# how many are detected only in humans?
hums <- psgndf %>%
  group_by(Genus, Species) %>%
  summarize(tot = sum(Abundance)) %>%
  ungroup() %>%
  mutate(pres = if_else(tot > 1, 1, 0)) %>%
  dplyr::select(-tot) %>%
  pivot_wider(names_from = Species, values_from = pres) %>%
  mutate(humanonly = if_else(homo_sapiens == 1 & cercocebus_agilis == 0 & gorilla_beringei == 0 & gorilla_gorilla == 0 & indri_indri == 0 & pan_troglodytes == 0 & papio_cynocephalus == 0 & procolobus == 0, "human", "no")) %>% filter(humanonly == "human") %>% pull(Genus)

# what is the relative abundance of these?
humdf <- psg.named %>% 
  tax_transform("compositional") %>% 
  tax_select(hums, strict_matches = TRUE) %>% 
  ps_filter(Species == "homo_sapiens") %>% 
  ps_melt() %>% 
  group_by(OTU) %>% 
  get_summary_stats(Abundance, type = "mean_sd")
summary(humdf$mean)

# what is the aggregated abundance of all of these fungi per human?
humtot <- psg.named %>% 
  #tax_transform("compositional") %>% 
  tax_select(hums, strict_matches = TRUE) %>% 
  ps_filter(Species == "homo_sapiens") %>%
  ps_melt() %>% 
  group_by(Sample) %>% 
  summarize(tot = sum(Abundance)) %>% 
  left_join(psg.named %>% 
              ps_mutate(seqdepth = sample_sums(psg.named)) %>% 
              samdat_tbl() %>% dplyr::select(Sample = .sample_name, seqdepth)) %>% 
  mutate(totabund = tot / seqdepth)
summary(humtot$totabund)
humtot %>% get_summary_stats(totabund)


# what is the prevalence of these fungi among humans?
humdfprev <- psg.named %>% 
  tax_transform("pa") %>% 
  tax_select(hums, strict_matches = TRUE) %>% 
  ps_filter(Species == "homo_sapiens") %>% 
  ps_melt() %>% 
  group_by(OTU) %>% 
  summarize(tot = sum(Abundance)) %>% 
  mutate(prev = tot / 33)
summary(humdfprev$prev) # mean 11% prev, min 6% and max 45%

# how many of these genera does each human have?
humdfcount <- psg.named %>% 
  tax_transform("pa") %>% 
  tax_select(hums, strict_matches = TRUE) %>% 
  ps_filter(Species == "homo_sapiens") %>% 
  ps_melt() %>% 
  group_by(Sample) %>% 
  summarize(tot = sum(Abundance)) %>% 
  mutate(totpresent = tot / 54)
summary(humdfcount$totpresent)
summary(humdfcount$tot)

## ---- chi squ test of unique genera ----

## get data
hum <- psgndf %>%
  group_by(Genus, Species) %>%
  summarize(tot = sum(Abundance)) %>%
  ungroup() %>%
  mutate(pres = if_else(tot > 1, 1, 0)) %>%
  dplyr::select(-tot)

# summarize: for each Species, how many Genus are unique? (not found in any other Species)
unique_genera <- hum %>%
  group_by(Genus) %>%
  summarize(total_species = sum(pres)) %>%
  left_join(hum, by = "Genus") %>%
  filter(total_species == 1 & pres == 1) %>%
  group_by(Species) %>%
  summarize(unique_genera = n()) %>%
  arrange(desc(unique_genera)) %>%
  # what is the percentage of total genera?
  mutate(perc = unique_genera / 515)

# compare the number of unique_genera between Species; significant?
fisher.test(unique_genera$Species, unique_genera$unique_genera, simulate.p.value = TRUE)
# turn this into a contingency table with the number of unique genera and the total number of genera (515)
ctab <- matrix(c(unique_genera$unique_genera, 515 - unique_genera$unique_genera), byrow = TRUE, nrow = 2)
rownames(ctab) <- c("Unique", "Non-Unique")
colnames(ctab) <- unique_genera$Species
chisq.test(ctab) # significant
pairwise.prop.test(t(ctab), p.adjust.method = "bonferroni")
# humans are significant to mangabey, papio, gorilla beringei, and procolobus


## ---- model each human-NHP comparison ----

# set up function
myModel <- function(phylo, taxa, comp_species) {
  taxtest <- taxatree_models2stats(phylo %>%
                                     ps_filter(Species %in% comp_species) %>%
                                     tax_select(taxa, strict_matches = TRUE) %>%
                                     tax_transform("compositional", rank = "unique") %>%
                                     # tax_filter(min_prevalence = 0.1, use_counts = TRUE) %>%
                                     taxatree_models(
                                       type = lm,
                                       trans = "log2", trans_args = list(zero_replace = "halfmin"),
                                       ranks = "unique", variables = "Species"
                                     )) %>%
    taxatree_stats_get() %>%
    mutate(comparisons = paste(comp_species, collapse = ":"))
  return(taxtest)
}


### human-v-pan troglodytes
pt_mod <- myModel(
  psg.named, prevg %>%
    filter(Species %in% c("homo_sapiens", "pan_troglodytes")) %>%
    filter(prevalence > 0.2) %>%
    distinct(OTU) %>%
    pull(OTU),
  c("homo_sapiens", "pan_troglodytes")
)

## human-v-gorilla gorilla
gg_mod <- myModel(
  psg.named, prevg %>%
    filter(Species %in% c("homo_sapiens", "gorilla_gorilla")) %>%
    filter(prevalence > 0.2) %>%
    distinct(OTU) %>%
    pull(OTU),
  c("homo_sapiens", "gorilla_gorilla")
)

## human-v-gorilla beringei
gb_mod <- myModel(
  psg.named, prevg %>%
    filter(Species %in% c("homo_sapiens", "gorilla_beringei")) %>%
    filter(prevalence > 0.2) %>%
    distinct(OTU) %>%
    pull(OTU),
  c("homo_sapiens", "gorilla_beringei")
)

## human-v-cercocebus agilis
ca_mod <- myModel(
  psg.named, prevg %>%
    filter(Species %in% c("homo_sapiens", "cercocebus_agilis")) %>%
    filter(prevalence > 0.2) %>%
    distinct(OTU) %>%
    pull(OTU),
  c("homo_sapiens", "cercocebus_agilis")
)

## human-v-indri indri
ii_mod <- myModel(
  psg.named, prevg %>%
    filter(Species %in% c("homo_sapiens", "indri_indri")) %>%
    filter(prevalence > 0.2) %>%
    distinct(OTU) %>%
    pull(OTU),
  c("homo_sapiens", "indri_indri")
)

## human-v-papio cynocephalus
pc_mod <- myModel(
  psg.named, prevg %>%
    filter(Species %in% c("homo_sapiens", "papio_cynocephalus")) %>%
    filter(prevalence > 0.2) %>%
    distinct(OTU) %>%
    pull(OTU),
  c("homo_sapiens", "papio_cynocephalus")
)

## human-v-procolobus
pr_mod <- myModel(
  psg.named, prevg %>%
    filter(Species %in% c("homo_sapiens", "procolobus")) %>%
    filter(prevalence > 0.2) %>%
    distinct(OTU) %>%
    pull(OTU),
  c("homo_sapiens", "procolobus")
)

#### add all together
all_mods <- rbind(pt_mod, gg_mod, gb_mod, ca_mod, ii_mod, pc_mod, pr_mod)
# FDR adjust
all_mods <- all_mods %>%
  mutate(padj = p.adjust(p.value, method = "fdr"), padj.bon = p.adjust(p.value, method = "bonferroni"))

# get sig stats
sig_mods <- all_mods %>%
  filter(padj < 0.05) # 84 total
length(unique(sig_mods$taxon)) # 173 total taxa are significant

# how many are differentiating between more than one NHP?
all_mods %>%
  filter(padj < 0.05) %>%
  group_by(taxon) %>%
  count() %>%
  arrange(desc(n))

# how many taxon were tested?
length(unique(all_mods$taxon))

# get this list of taxa that were different between all 7 NHP
sig_taxa <- all_mods %>%
  filter(padj < 0.05) %>%
  group_by(taxon) %>%
  count() %>%
  filter(n >= 6) %>%
  pull(taxon)

# wrangle output
all_mods <- all_mods %>%
  mutate(comp = case_when(
    comparisons %in% "homo_sapiens:pan_troglodytes" ~ paste(taxon, "P_troglodytes"),
    comparisons %in% "homo_sapiens:gorilla_gorilla" ~ paste(taxon, "G_gorilla"),
    comparisons %in% "homo_sapiens:gorilla_beringei" ~ paste(taxon, "G_beringei"),
    comparisons %in% "homo_sapiens:cercocebus_agilis" ~ paste(taxon, "C_agilis"),
    comparisons %in% "homo_sapiens:indri_indri" ~ paste(taxon, "I_indri"),
    comparisons %in% "homo_sapiens:papio_cynocephalus" ~ paste(taxon, "P_cynocephalus"),
    comparisons %in% "homo_sapiens:procolobus" ~ paste(taxon, "P_gordonum")
  ))


# calculate log2 change
lg <- psg.named %>%
  tax_select(sig_taxa, strict_matches = TRUE) %>%
  tax_transform("compositional") %>%
  tax_transform("log2", zero_replace = "halfmin", chain = TRUE) %>%
  ps_melt() 

# get the means
lgm <- lg %>%
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

lgm <- lgm %>%
  # order fungal genera
  mutate(gen.order = factor(OTU, ordered = TRUE, levels = lgm %>% arrange(desc(fc)) %>% distinct(OTU) %>% pull(OTU)))

# write out for summary statistics
sums <- lgm %>%
  dplyr::select(SciName, Genus = OTU, fc, estimate, std.error, statistic, padj) %>%
  arrange(SciName, Genus)
# write out
write.table(sums, file = "data/results_diffrelabund_humanvNHP.txt", sep = "\t", row.names = FALSE)




