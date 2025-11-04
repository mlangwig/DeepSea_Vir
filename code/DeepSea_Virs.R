library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(stringr)
library(grid)

############################### Input ########################################

# Read img input
imgvr4 <- read.delim("/scratch/langwig/Projects/DeepSea_VirusReview/IMGVR4_all_Sequence_information_Marine.tsv", header = FALSE)

# Assign column names
colnames(imgvr4) <- c("UVIG", "Taxon_oid", "Scaffold_oid", "Coordinates ('whole' if the UViG is the entire contig)", "Ecosystem classification", "vOTU", "Length", "Topology", "geNomad score", "Confidence", "Estimated completeness", "Estimated contamination", "MIUViG quality", "Gene content (total genes;cds;tRNA;geNomad marker)", "Taxonomic classification", "Taxonomic classification method", "Host taxonomy prediction", "Host prediction method", "Sequence origin (doi)")

############################### Processing ########################################

# Remove non-deep sea data
imgvr4 <- imgvr4 %>%
library(dplyr)

imgvr4_filt <- imgvr4 %>%
    filter(!grepl("Coastal|Continental margin|Deep subsurface|Epipelagic|Fjord|Fossil|Harbor|Inlet|Intertidal zone|Marine lake|Mud volcano|Neritic zone|River plume|Sea ice|Seaweed|Subtidal zone|Unclassified|Volcanic|Wetlands", `Ecosystem classification`))
  #filter(!`Ecosystem classification` %in% c("Coastal", "Continental margin", "Deep subsurface", "Epipelagic", "Fjord", "Fossil", "Harbor", "Inlet", "Intertidal zone", "Marine lake", "Mud volcano", "Neritic zone", "Neritic zone/Coastal water", "River plume", "Sea ice", "Seaweed", "Subtidal zone", "Unclassified", "Volcanic", "Wetlands"))
    #filter(!grepl("Coastal|Continental margin|Deep subsurface|Epipelagic|Fjord|Fossil|Harbor|Inlet|Intertidal|Neritic|River plume|Sea ice|Seaweed|Wetlands|Volcanic|Unclassified|Subtidal zone|", `Ecosystem classification`))