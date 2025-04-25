library(tidyverse)
library(phyloseq)
library(microViz)


# read in spreadsheet
# This can be replaced with Supplementary Table 1 from paper
meta <- readxl::read_xlsx(
  path = "Data_PrimateITS_v1.xlsx",
  sheet = "Species_List"
)

# read in each runids file and make master list of Run accessions
runids <- list.files("all_primate", pattern = "runids_*", full.names = TRUE)
rundf <- data.frame()
for (i in 1:length(runids)) {
  subdf <- read.table(runids[i], header = FALSE)
  subdf <- subdf %>%
    mutate(Accession = str_extract(runids[i], "PRJ\\w+(_\\w+)?")) %>%
    dplyr::rename(Run = V1) %>%
    mutate(is.subset = if_else(str_detect(runids[i], "subset"), TRUE, FALSE))
  rundf <- bind_rows(rundf, subdf)
}

# join the metadata for accessions that are single-species datasets (only Indri)
run39443 <- meta %>%
  filter(Accession == "PRJEB39443") %>%
  left_join(rundf, by = "Accession")

# get metadata for accessions with multiple species
run37770 <- read.table("all_primate/PRJEB37770_Zenodo/metadata_PRJEB37770_Zenodo.txt", sep = "\t", header = TRUE) %>%
  # keep only columns of interest
  dplyr::select(Run, Species) %>%
  mutate(Run = str_remove(Run, "ITS3_ITS4-")) %>%
  mutate(Run = str_remove(Run, "_[A-Z]+$")) %>%
  left_join(rundf, by = "Run") %>% # 158 samples and 2 host species; matches
  mutate(
    Species = tolower(Species),
    Accession = "PRJEB37770"
  ) %>%
  left_join(meta)

## this was provided by the study authors
run686661 <- readxl::read_xlsx("PRJNA686661_metadata.xlsx") %>%
  # remove unknown samples
  filter(!FileID == "Unknown") %>%
  # rename(Run = FileID) %>%
  inner_join(rundf %>% mutate(Run = str_remove(Run, "_S(\\d){1,4}$")), by = c("FileID" = "Run"), keep = TRUE) %>%
  # wrangle
  mutate(Species = case_when(
    Group %in% "BaAka-Human" ~ "homo sapiens",
    Group %in% "Bantu-Human" ~ "homo sapiens",
    Group %in% "Chimps" ~ "pan troglodytes schweinfurthii",
    Group %in% "Mangabey" ~ "cercocebus agilis",
    Group %in% "Mountain Gorilla" ~ "gorilla beringei",
    Group %in% "Western Lowland Gorilla" ~ "gorilla gorilla",
  )) %>%
  # keep only columns of interest
  dplyr::select(Run = FileID, Species, Accession) %>%
  left_join(meta)

#### join all metadata
allmeta <- bind_rows(run39443, run37770, run686661) %>%
  drop_na(Run)

### ----- add study metadata -----


# read in study-specific metadata
# can be replaced with Supplementary Table 2
studdf <- readxl::read_xlsx(path = "all_primate/Data_PrimateITS_v1.xlsx", sheet = "Samples_List") %>%
  # do some wrangling
  filter(SRA_Accession %in% allmeta$Accession)

# what different primer sets were used?
studdf <- studdf %>% mutate(primerset = paste(FWD_Primer_Name, REV_Primer_Name, sep = ":"))
studdf %>%
  distinct(SRA_Accession, `Sequencing type`, primerset) %>%
  group_by(primerset) %>%
  count()

# get primer info to use
pdf <- studdf %>%
  dplyr::select(Accession = SRA_Accession, primerset, Sequencing_Type = `Sequencing type`) %>%
  distinct() %>%
  mutate(Sequencing_Type = if_else(str_detect(Sequencing_Type, "18S"), "18S", Sequencing_Type))
allmeta <- allmeta %>%
  left_join(pdf)

### ---- write out ----


# save
write.table(allmeta, file = "data/all_primate_metadata.txt", sep = "\t", row.names = FALSE)
# save just runIDS for bioinformatics
write.table(allmeta$Run,
            file = "/data/all_primate_runids.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t"
)
