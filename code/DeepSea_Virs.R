library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(stringr)
library(grid)
library(readxl)

############################### Input ########################################

# Read img input
imgvr4 <- read.delim("/scratch/langwig/Projects/DeepSea_VirusReview/IMGVR4_all_Sequence_information_Marine.tsv", header = FALSE)

# Assign column names
colnames(imgvr4) <- c("UVIG", "Taxon_oid", "Scaffold_oid", "Coordinates ('whole' if the UViG is the entire contig)", "Ecosystem classification", "vOTU", "Length", "Topology", "geNomad score", "Confidence", "Estimated completeness", "Estimated contamination", "MIUViG quality", "Gene content (total genes;cds;tRNA;geNomad marker)", "Taxonomic classification", "Taxonomic classification method", "Host taxonomy prediction", "Host prediction method", "Sequence origin (doi)")

# GOLD DB mapping files - excel sheet downloaded from https://gold.jgi.doe.gov/downloads
gold_taxid <- read_excel("/scratch/langwig/Projects/DeepSea_VirusReview/goldData.xlsx", sheet = "Analysis Project")
gold_metadata <- read_excel("/scratch/langwig/Projects/DeepSea_VirusReview/goldData.xlsx", sheet = "Sequencing Project") #615,448 rows 

############################### Processing ########################################

# Remove non-deep sea data
imgvr4_filt <- imgvr4 %>%
    filter(!grepl("Coastal|Continental margin|Deep subsurface|Epipelagic|Fjord|Fossil|Harbor|Inlet|Intertidal zone|Marine lake|Mud volcano|Neritic zone|River plume|Sea ice|Seaweed|Subtidal zone|Unclassified|Volcanic|Wetlands", `Ecosystem classification`))
  #filter(!`Ecosystem classification` %in% c("Coastal", "Continental margin", "Deep subsurface", "Epipelagic", "Fjord", "Fossil", "Harbor", "Inlet", "Intertidal zone", "Marine lake", "Mud volcano", "Neritic zone", "Neritic zone/Coastal water", "River plume", "Sea ice", "Seaweed", "Subtidal zone", "Unclassified", "Volcanic", "Wetlands"))
    #filter(!grepl("Coastal|Continental margin|Deep subsurface|Epipelagic|Fjord|Fossil|Harbor|Inlet|Intertidal|Neritic|River plume|Sea ice|Seaweed|Wetlands|Volcanic|Unclassified|Subtidal zone|", `Ecosystem classification`))

# Second pass after some manual checking of project IDs
# imgvr4_filt <- imgvr4 %>%
#     filter(!grepl("Strait|Gulf", `Ecosystem classification`))

############################### Processing to get more specific than Oceanic ########################################

# Subset just oceanic - 1,476,085 rows
imgvr4_oceanic <- imgvr4_filt %>%
    filter(`Ecosystem classification` == "Environmental;Aquatic;Marine;Oceanic")

# Make tax id col character rather than string
imgvr4_oceanic$Taxon_oid <- as.character(imgvr4_oceanic$Taxon_oid)

# Match taxid to table with GOLD ID
imgvr4_oceanic <- gold_taxid %>%
  dplyr::select(c("AP IMG TAXON ID", "AP PROJECT GOLD IDS")) %>%
  right_join(imgvr4_oceanic, by = c("AP IMG TAXON ID" = "Taxon_oid"))

# See how successful that was
sum(is.na(imgvr4_oceanic$`AP PROJECT GOLD IDS`))
#1 no match?

# check gold metadata
gold_metadata_filt <- gold_metadata %>%
  dplyr::select(c("PROJECT GOLD ID", "PROJECT NAME")) %>%
  unique()

# Match B and R cols in excel
imgvr4_oceanic_tst <- gold_metadata_filt %>%
  right_join(imgvr4_oceanic, by = c("PROJECT GOLD ID" = "AP PROJECT GOLD IDS"))

# See how successful that was
sum(is.na(imgvr4_oceanic_tst$`PROJECT NAME`))
#34,981 with no project name

# checking number of occurrences
sum(str_detect(imgvr4_oceanic_tst$`PROJECT NAME`, "Deepwater Horizon"), na.rm = TRUE)

envmts <- as.data.frame(unique(imgvr4_oceanic_tst$`PROJECT NAME`))

########### Deciding who to remove and rename ########################



############################### Create plotting table ########################################

# Group and sum for counts
plot <- imgvr4_filt %>%
  group_by(`Ecosystem classification`) %>%
  summarise(UVIG_count = n_distinct(UVIG))
# Original sum here I have count 5 for Environmental;Aquatic;Marine; and count 2 for Environmental;Aquatic;Marine;Marine basin floor
# Note that I am arbitrarily combining these with Environmental;Aquatic;Marine;Oceanic and Environmental;Aquatic;Marine;Benthic respectively

# Add 5 and 2 to get rid of small categories that are indistinguishable
plot <- plot %>%
  mutate(`Ecosystem classification` = if_else(
    `Ecosystem classification` == "Environmental;Aquatic;Marine;",
    "Environmental;Aquatic;Marine;Oceanic",
    `Ecosystem classification`)) %>%
  group_by(`Ecosystem classification`) %>%
  summarise(UVIG_count = sum(UVIG_count), .groups = "drop")

plot <- plot %>%
  mutate(`Ecosystem classification` = if_else(
    `Ecosystem classification` == "Environmental;Aquatic;Marine;Marine basin floor",
    "Environmental;Aquatic;Marine;Benthic",
    `Ecosystem classification`)) %>%
  group_by(`Ecosystem classification`) %>%
  summarise(UVIG_count = sum(UVIG_count), .groups = "drop")

# Get the most specific name
plot <- plot %>%
  mutate(Specific_env = str_replace(`Ecosystem classification`, ";$", "")) %>%  # remove trailing ;
  mutate(Specific_env = word(Specific_env, -1, sep = fixed(";"))) %>%           # take the last level
  arrange(UVIG_count)        

################################ Lollipop plot #########################################

ggplot(plot, aes(x = UVIG_count, y = reorder(Specific_env, UVIG_count))) +
  geom_segment(aes(x = 0, xend = UVIG_count, yend = Specific_env), color = "grey60") +
  geom_point(size = 3, color = "steelblue") +
  labs(
    x = "UVIG Count",
    y = "Ecosystem",
    title = ""
  ) +
  theme_minimal(base_size = 13) +
  scale_x_continuous(expand = c(0,0), labels = comma)
