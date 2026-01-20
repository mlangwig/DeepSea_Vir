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

# IMG/VR custom input
imgvr4 <- read.delim("/scratch/langwig/Projects/DeepSea_VirusReview/DeepSea_Vir/data/exported_img_data.tsv", header = TRUE)

# # Read img input
# imgvr4 <- read.delim("/scratch/langwig/Projects/DeepSea_VirusReview/IMGVR4_all_Sequence_information_Marine.tsv", header = FALSE)

# # Assign column names
# colnames(imgvr4) <- c("UVIG", "Taxon_oid", "Scaffold_oid", "Coordinates ('whole' if the UViG is the entire contig)", "Ecosystem classification", "vOTU", "Length", "Topology", "geNomad score", "Confidence", "Estimated completeness", "Estimated contamination", "MIUViG quality", "Gene content (total genes;cds;tRNA;geNomad marker)", "Taxonomic classification", "Taxonomic classification method", "Host taxonomy prediction", "Host prediction method", "Sequence origin (doi)")

# # GOLD DB mapping files - excel sheet downloaded from https://gold.jgi.doe.gov/downloads
# gold_taxid <- read_excel("/scratch/langwig/Projects/DeepSea_VirusReview/goldData.xlsx", sheet = "Analysis Project")
# gold_metadata <- read_excel("/scratch/langwig/Projects/DeepSea_VirusReview/goldData.xlsx", sheet = "Sequencing Project") #615,448 rows 

############################### Processing ########################################

# Remove rows where no viruses ID'd
imgvr4_filt <- imgvr4 %>%
  filter(`Predicted.Viruses` > 0)

# how many viruses >200 m depth?
sum(imgvr4_filt$Predicted.Viruses)
#757,054

# how many have no depth listed? ## DECIDE IF YOU JUST WANT TO REMOVE THESE VIRUSES? ##
sum(is.na(imgvr4_filt$Depth.In.Meters) | imgvr4_filt$Depth.In.Meters == "")
# 65 projects
# how many viral genomes does this encompass?
imgvr4_filt %>%
  filter(is.na(Depth.In.Meters) | Depth.In.Meters == "") %>%
  summarise(total = sum(Predicted.Viruses, na.rm = TRUE))
# 191 viruses - negligible

# # OMZs
# tst_omz <- imgvr4_filt %>%
#   mutate(Specific.Ecosystem = if_else(
#       grepl("OMZ|oxygen minimum zones|oxygen minimum zone|Oxygen Minimum Zone|Oxygen Minumum", Genome.Name...Sample.Name),
#       "Hypoxic zone", Specific.Ecosystem))

# tst_omz <- tst_omz %>%
#   filter(Specific.Ecosystem == "Hypoxic zone")

# sum(tst_omz$Predicted.Viruses)

############################### Unclassified ########################################

# filter unclassified from img table, reclassify according to sample name
imgvr4_filt_unclass <- imgvr4_filt %>%
  filter(`Ecosystem.Subtype` == "Unclassified")

# change depth to numeric
imgvr4_filt_unclass$Depth.In.Meters <- as.numeric(imgvr4_filt_unclass$Depth.In.Meters)

# Tara Oceans and seawater
# change to more specific by depth
imgvr4_filt_unclass <- imgvr4_filt_unclass %>%
  mutate(
    # 1 Tara Oceans samples as Pelagic
    Ecosystem.Subtype = if_else(
      grepl("Tara Oceans|Seawater|seawater|pool|Red Sea|Sea|Ocean", Genome.Name...Sample.Name),
      "Pelagic", Ecosystem.Subtype
    )
  ) %>%
  mutate(
    # 2 change Pelagic by depth
    Ecosystem.Subtype = case_when(
      Ecosystem.Subtype == "Pelagic" & 
        Depth.In.Meters >= 200 & Depth.In.Meters < 1000 ~ "Mesopelagic",

      Ecosystem.Subtype == "Pelagic" & 
        Depth.In.Meters >= 1000 & Depth.In.Meters < 4000 ~ "Bathypelagic",

      Ecosystem.Subtype == "Pelagic" & 
        Depth.In.Meters >= 4000 & Depth.In.Meters <= 6000 ~ "Abyssopelagic",

      Ecosystem.Subtype == "Pelagic" & 
        Depth.In.Meters > 6000 ~ "Hadopelagic",

      TRUE ~ Ecosystem.Subtype
    )
  )

# Microbial mats
imgvr4_filt_unclass <- imgvr4_filt_unclass %>%
  mutate(Ecosystem.Subtype = if_else(
      grepl("Microbial mat", Genome.Name...Sample.Name),
      "Microbial mats", Ecosystem.Subtype))

# OMZs
imgvr4_filt_unclass <- imgvr4_filt_unclass %>%
  mutate(Ecosystem.Subtype = if_else(
      grepl("OMZs|low-oxygen", Study.Name),
      "Hypoxic zone", Ecosystem.Subtype))

imgvr4_filt_unclass <- imgvr4_filt_unclass %>%
  mutate(Specific.Ecosystem = if_else(
      grepl("OMZ|oxygen minimum zones|oxygen minimum zone|Oxygen Minimum Zone|Oxygen Minumum", Genome.Name...Sample.Name),
      "Hypoxic zone", Specific.Ecosystem))

# Microbial isolates
imgvr4_filt_unclass <- imgvr4_filt_unclass %>%
  mutate(Ecosystem.Subtype = if_else(
      grepl("Genome sequencing|Genome Sequencing|genome sequencing|Strains|strains|bacteria sequencing|bacterium|Shewanella", Study.Name),
      "Microbial isolate", Ecosystem.Subtype))

imgvr4_filt_unclass <- imgvr4_filt_unclass %>%
  mutate(Ecosystem.Subtype = if_else(
      grepl("Candidatus|candidate|DSM|PV-1", Genome.Name...Sample.Name),
      "Microbial isolate", Ecosystem.Subtype))

imgvr4_filt_unclass <- imgvr4_filt_unclass %>%
  mutate(
    Ecosystem.Subtype = if_else(
      grepl("(?<![A-Za-z0-9])sp\\.(?![A-Za-z0-9])",
            Genome.Name...Sample.Name,
            ignore.case = TRUE,
            perl = TRUE),
      "Microbial isolate",
      Ecosystem.Subtype
    )
  )

# Sediment
imgvr4_filt_unclass <- imgvr4_filt_unclass %>%
  mutate(Ecosystem.Subtype = if_else(
      grepl("sediment", Study.Name),
      "Sediment", Ecosystem.Subtype))

# Group and sum for counts
plot_unclass <- imgvr4_filt_unclass %>%
  group_by(`Ecosystem.Subtype`) %>%
  summarise(UVIG_count = sum(Predicted.Viruses))

# ############################### Oceanic ########################################

# Working with Specific.Ecosystem col here instead of Ecosystem.Subtype like Unclassified

# Remove Unclassifieds from og
imgvr4_filt <- imgvr4_filt %>%
  filter(!`Ecosystem.Subtype` == "Unclassified")

# Parse the Oceanic into more specifics - IF oceanic THEN use specific.ecosystem?
imgvr4_filt_oceanic <- imgvr4_filt %>%
  filter(`Ecosystem.Subtype` == "Oceanic")

# change depth to numeric
imgvr4_filt_oceanic$Depth.In.Meters <- as.numeric(imgvr4_filt_oceanic$Depth.In.Meters)

# Water column
# change to more specific by depth
imgvr4_filt_oceanic <- imgvr4_filt_oceanic %>%
  mutate(
    # 1 change water columns to pelagic
    Specific.Ecosystem = if_else(
      grepl("Malaspina|Tara Oceans|Seawater|seawater|pool|Red Sea|Sea|Ocean|Gyre|Pond", Genome.Name...Sample.Name),
      "Pelagic", Specific.Ecosystem
    )
  ) %>%
  mutate(
    # 2 change Pelagic by depth
    Specific.Ecosystem = case_when(
      Specific.Ecosystem == "Pelagic" & 
        Depth.In.Meters >= 200 & Depth.In.Meters < 1000 ~ "Mesopelagic",

      Specific.Ecosystem == "Pelagic" & 
        Depth.In.Meters >= 1000 & Depth.In.Meters < 4000 ~ "Bathypelagic",

      Specific.Ecosystem == "Pelagic" & 
        Depth.In.Meters >= 4000 & Depth.In.Meters <= 6000 ~ "Abyssopelagic",

      Specific.Ecosystem == "Pelagic" & 
        Depth.In.Meters > 6000 ~ "Hadopelagic",

      TRUE ~ Specific.Ecosystem
    )
  )

# Microbial isolates
imgvr4_filt_oceanic <- imgvr4_filt_oceanic %>%
  mutate(Specific.Ecosystem = if_else(
      grepl("isolated|Strains|strains|bacteria sequencing|bacterium", Study.Name),
      "Microbial isolate", Specific.Ecosystem))

imgvr4_filt_oceanic <- imgvr4_filt_oceanic %>%
  mutate(Specific.Ecosystem = if_else(
      grepl("Candidatus|candidate|DSM|PV-1", Genome.Name...Sample.Name),
      "Microbial isolate", Specific.Ecosystem))

imgvr4_filt_oceanic <- imgvr4_filt_oceanic %>%
  mutate(
    Specific.Ecosystem = if_else(
      grepl("(?<![A-Za-z0-9])sp\\.(?![A-Za-z0-9])",
            Genome.Name...Sample.Name,
            ignore.case = TRUE,
            perl = TRUE),
      "Microbial isolate",
      Specific.Ecosystem
    )
  )


# OMZs
imgvr4_filt_oceanic <- imgvr4_filt_oceanic %>%
  mutate(Specific.Ecosystem = if_else(
      grepl("OMZ|oxygen minimum zones|oxygen minimum zone|Oxygen Minimum Zone|Oxygen Minumum", Genome.Name...Sample.Name),
      "Hypoxic zone", Specific.Ecosystem))

# Hydro vent
imgvr4_filt_oceanic <- imgvr4_filt_oceanic %>%
  mutate(Specific.Ecosystem = if_else(
      grepl("vent plume", Study.Name),
      "Hydrothermal vents", Specific.Ecosystem))



# Group and sum for counts
plot_oceanic <- imgvr4_filt_oceanic %>%
  group_by(`Specific.Ecosystem`) %>%
  summarise(UVIG_count = sum(Predicted.Viruses))



# Classify the NCBI.Domain archaea and bacteria as deep-sea isolates?


############################### Create plotting table ########################################

# Remove Unclassifieds from og
imgvr4_filt <- imgvr4_filt %>%
  filter(!`Ecosystem.Subtype` == "Unclassified")

# Remove Oceanics from og
imgvr4_filt <- imgvr4_filt %>%
  filter(!`Ecosystem.Subtype` == "Oceanic")

# Group and sum for counts
plot <- imgvr4_filt %>%
  group_by(`Ecosystem.Subtype`) %>%
  summarise(UVIG_count = sum(Predicted.Viruses))

# rename oceanic plot col name so can combine
plot_oceanic <- plot_oceanic %>%
  rename(`Ecosystem.Subtype` = `Specific.Ecosystem`) # new = old

# combine og, unclassifieds, and oceanics
plot_final <- rbind(plot, plot_unclass, plot_oceanic)

# re sum
plot_final <- plot_final %>%
  group_by(`Ecosystem.Subtype`) %>%
  summarise(UVIG_count = sum(UVIG_count))

# remove low counts
plot_final <- plot_final %>%
  filter(UVIG_count > 2000)

lolly <- ggplot(plot_final, aes(x = UVIG_count, y = reorder(Ecosystem.Subtype, UVIG_count))) +
  geom_segment(aes(x = 0, xend = UVIG_count, yend = Ecosystem.Subtype), color = "#000000", linewidth = 0.7) +
  geom_point(size = 5, color = "#209bcf") +
  # geom_point(
  #   data = subset(plot, Cultured == "no"),
  #   aes(x = UVIG_count + 0.05 * max(UVIG_count)),
  #   shape = 17, size = 4, color = "red"
  # ) +
  labs(
    x = "UVIG Count",
    y = "Ecosystem",
    title = ""
  ) +
  theme_minimal(base_size = 18) +
  scale_x_continuous(expand = c(.01,0), labels = comma) +
  theme(
    panel.grid.major.y = element_blank(),  # remove horizontal gridlines
    panel.grid.minor.y = element_blank(),
    
    panel.grid.major.x = element_line(linetype = "dashed"),  # dashed vertical gridlines
    panel.grid.minor.x = element_blank()
  )
lolly

ggsave("/scratch/langwig/Projects/DeepSea_VirusReview/DeepSea_Vir/visuals/UVIGs_DeepSea_IMGVR.png", 
lolly, dpi = 500, bg = "transparent") #change height to 15 when mod'ing legend

############################### Old Processing ########################################

# # # See what the Strait environment is about
# # imgvr4_filt_Strait <- imgvr4 %>%
# #   filter(`Ecosystem classification` == "Environmental;Aquatic;Marine;Strait")

# # # See what the Strait environment is about
# # imgvr4_filt_Aquifer <- imgvr4 %>%
# #   filter(`Ecosystem classification` == "Environmental;Aquatic;Marine;Aquifer")

# # # See what the Hypoxic zone environment is about
# # imgvr4_filt_Hypoxic <- imgvr4 %>%
# #   filter(`Ecosystem classification` == "Environmental;Aquatic;Marine;Hypoxic zone")

# # Remove non-deep sea data
# imgvr4_filt <- imgvr4 %>%
#     filter(!grepl("Coastal|Strait|Aquifer|aquifer|Continental margin|Deep subsurface|Epipelagic|Fjord|Fossil|Harbor|Inlet|Intertidal zone|Marine lake|Mud volcano|Neritic zone|River plume|Sea ice|Seaweed|Subtidal zone|Unclassified|Volcanic|Wetlands", `Ecosystem classification`))
#   #filter(!`Ecosystem classification` %in% c("Coastal", "Continental margin", "Deep subsurface", "Epipelagic", "Fjord", "Fossil", "Harbor", "Inlet", "Intertidal zone", "Marine lake", "Mud volcano", "Neritic zone", "Neritic zone/Coastal water", "River plume", "Sea ice", "Seaweed", "Subtidal zone", "Unclassified", "Volcanic", "Wetlands"))
#     #filter(!grepl("Coastal|Continental margin|Deep subsurface|Epipelagic|Fjord|Fossil|Harbor|Inlet|Intertidal|Neritic|River plume|Sea ice|Seaweed|Wetlands|Volcanic|Unclassified|Subtidal zone|", `Ecosystem classification`))

# # Second pass after some manual checking of project IDs
# # imgvr4_filt <- imgvr4 %>%
# #     filter(!grepl("Strait|Gulf", `Ecosystem classification`))

# ############################### Processing to get more specific than Oceanic ########################################

# # Subset just oceanic - 1,476,085 rows
# imgvr4_oceanic <- imgvr4_filt %>%
#     filter(`Ecosystem classification` == "Environmental;Aquatic;Marine;Oceanic")

# # Make tax id col character rather than numeric
# imgvr4_oceanic$Taxon_oid <- as.character(imgvr4_oceanic$Taxon_oid)

# # Match taxid to table with GOLD ID
# imgvr4_oceanic <- gold_taxid %>%
#   dplyr::select(c("AP IMG TAXON ID", "AP PROJECT GOLD IDS")) %>%
#   right_join(imgvr4_oceanic, by = c("AP IMG TAXON ID" = "Taxon_oid"))

# # See how successful that was
# sum(is.na(imgvr4_oceanic$`AP PROJECT GOLD IDS`))
# #1 no match?

# # check gold metadata
# gold_metadata_filt <- gold_metadata %>%
#   dplyr::select(c("PROJECT GOLD ID", "PROJECT NAME")) %>%
#   unique()

# # Match B and R cols in excel
# imgvr4_oceanic_tst <- gold_metadata_filt %>%
#   right_join(imgvr4_oceanic, by = c("PROJECT GOLD ID" = "AP PROJECT GOLD IDS"))

# # See how successful that was
# sum(is.na(imgvr4_oceanic_tst$`PROJECT NAME`))
# #34,981 with no project name

# # How many times do project names occur
# project_counts <- imgvr4_oceanic_tst %>%
#   count(`PROJECT NAME`, name = "count") %>%
#   arrange(desc(count))
# # This showed me that a lot of the top counts for Oceanic are what I think is Tara Oceans and most of these are shallow

# # checking number of occurrences
# sum(str_detect(imgvr4_oceanic_tst$`PROJECT NAME`, "Deepwater Horizon"), na.rm = TRUE)
# sum(str_detect(imgvr4_oceanic_tst$`PROJECT NAME`, "oxygen minimum zone"), na.rm = TRUE)
# sum(str_detect(imgvr4_oceanic_tst$`PROJECT NAME`, regex("sediment", ignore_case = TRUE)), na.rm = TRUE)
# sum(str_detect(imgvr4_oceanic_tst$`PROJECT NAME`, "Seawater"), na.rm = TRUE)
# sum(str_detect(imgvr4_oceanic_tst$`PROJECT NAME`, regex("700m", ignore_case = TRUE)), na.rm = TRUE)
# sum(str_detect(imgvr4_oceanic_tst$`PROJECT NAME`, regex("JdFR", ignore_case = TRUE)), na.rm = TRUE)
# sum(str_detect(imgvr4_oceanic_tst$`PROJECT NAME`, regex(" 10m", ignore_case = TRUE)), na.rm = TRUE)

# envmts <- as.data.frame(unique(imgvr4_oceanic_tst$`PROJECT NAME`))
# colnames(envmts) <- "unique_env"

# ########### Deciding who to remove and rename ########################

# # clear cases and patterns to match for removal
# envmts_tst <- envmts %>%
#   filter(!str_detect(unique_env, regex("Prochlorococcus|phytoplankton|surface|plastic|mangrove|subsurface|Inlet|Oct19|AMALJGI-DNA-1|AMALJGI-DNA-2|AMALJGI-DNA-9|AMALJGI-DNA-10|_10m|5mSIP|Malaspina viral metaG Antarct|Costa Rica Dome|ArcticOcean|SST-1|SST-2|Lowphox|Arctic Ocean - MOSAiC| 10m|P04_10|SI03|_150m|_50m|CsCl metaG|[0-9]_ETSP_OMZ_AT|[0-9]B_ETSP_OMZ_AT", ignore_case = TRUE)))

# # clear cases and patterns to match for removal
# imgvr4_oceanic_metadata_filt <- imgvr4_oceanic_tst %>%
#   filter(!str_detect(`PROJECT NAME`, regex("Prochlorococcus|phytoplankton|surface|plastic|mangrove|subsurface|Inlet|Oct19|AMALJGI-DNA-1|AMALJGI-DNA-2|AMALJGI-DNA-9|AMALJGI-DNA-10|_10m|5mSIP|Malaspina viral metaG Antarct|Costa Rica Dome|ArcticOcean|SST-1|SST-2|Lowphox|Arctic Ocean - MOSAiC| 10m|P04_10|SI03|_150m|_50m|CsCl metaG|[0-9]_ETSP_OMZ_AT|[0-9]B_ETSP_OMZ_AT", ignore_case = TRUE)))

# # Manual searches of info above. Descriptions below...
# # Samples with names "Seawater viral communities from ... Oct19" are from 5 m depth
# # North Sea says 0 m
# # AMALJGI-DNA-9 and 10 are <200 m
# # JdFR seawater samples are all 2,000-3,000 m
# # Gs0067852 (Tara Oceans?) had very non specific env type and no way to drop the shallows easily so left as oceanic for removal
# # Costa Rica Dome are surface samples
# # ArcticOcean are all from Gs0118432, shallow
# # SST-1 and 2 from Gs0135243, surface samples
# # Lowphox all from Gs0067852, shallow
# # Arctic Ocean - MOSAiC from Gs0153906, shallow
# # [0-9]B_ETSP_OMZ_AT associated with shallow OMZ

# # # check logic
# # envmts_tst <- envmts %>%
# #   mutate(`Ecosystem classification` = case_when(
# #     str_detect(envmts$unique_env, regex("oxygen minimum zone", ignore_case = TRUE)) ~ "OMZ",
# #     str_detect(envmts$unique_env, regex("OMZ", ignore_case = TRUE)) ~ "OMZ",
# #     str_detect(envmts$unique_env, regex("oil-polluted", ignore_case = TRUE)) ~ "Oil polluted",
# #     str_detect(envmts$unique_env, regex("Deepwater Horizon", ignore_case = TRUE)) ~ "Oil polluted",
# #     str_detect(envmts$unique_env, regex("oil contamination", ignore_case = TRUE)) ~ "Oil polluted",
# #     str_detect(envmts$unique_env, regex("seawater", ignore_case = TRUE)) ~ "Pelagic zone",
# #     str_detect(envmts$unique_env, regex("Gulf", ignore_case = TRUE)) ~ "Gulf",
# #     str_detect(envmts$unique_env, regex("sediment", ignore_case = TRUE)) ~ "Sediment",
# #     str_detect(envmts$unique_env, regex("Sea bed", ignore_case = TRUE)) ~ "Benthic",
# #     #str_detect(envmts$unique_env, regex("aquifer", ignore_case = TRUE)) ~ "Aquifer",
# #     str_detect(envmts$unique_env, regex("AMALGI", ignore_case = TRUE)) ~ "Pelagic zone",
# #     str_detect(envmts$unique_env, regex("P26_", ignore_case = TRUE)) ~ "Pelagic zone",
# #     str_detect(envmts$unique_env, regex("P4_", ignore_case = TRUE)) ~ "Pelagic zone",
# #     str_detect(envmts$unique_env, regex("JdFR", ignore_case = TRUE)) ~ "Bathypelagic",
# #     str_detect(envmts$unique_env, regex("700m", ignore_case = TRUE)) ~ "Mesopelagic",
# #     str_detect(envmts$unique_env, regex("1000", ignore_case = TRUE)) ~ "Bathypelagic",
# #     str_detect(envmts$unique_env, regex("Malaspina", ignore_case = TRUE)) ~ "Bathypelagic",
# #     str_detect(envmts$unique_env, regex("MSP", ignore_case = TRUE)) ~ "Bathypelagic",
# #     str_detect(envmts$unique_env, regex("- MP", ignore_case = TRUE)) ~ "Bathypelagic",
# #     str_detect(envmts$unique_env, regex("500_MG", ignore_case = TRUE)) ~ "Mesopelagic",
# #     str_detect(envmts$unique_env, regex("200_MG", ignore_case = TRUE)) ~ "Mesopelagic",
# #     str_detect(envmts$unique_env, regex("WHOI_OMZ", ignore_case = TRUE)) ~ "Pelagic zone",
# #     str_detect(envmts$unique_env, regex("- LP-", ignore_case = TRUE)) ~ "Pelagic zone",
# #     TRUE ~ unique_env  # leave everything else as-is
# # ))

# # Then remove rows where only label Oceanic remains - these are ones I know to be shallow but was no easy way to match them all without omitting samples from same project that were deeper
# ####
# imgvr4_oceanic_metadata_filt_mod <- imgvr4_oceanic_metadata_filt %>%
#   mutate(`Ecosystem classification` = case_when(
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("oxygen minimum zone", ignore_case = TRUE)) ~ "OMZ",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("OMZ", ignore_case = TRUE)) ~ "OMZ",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("oil-polluted", ignore_case = TRUE)) ~ "Oil polluted",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("Deepwater Horizon", ignore_case = TRUE)) ~ "Oil polluted",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("oil contamination", ignore_case = TRUE)) ~ "Oil polluted",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("seawater", ignore_case = TRUE)) ~ "Pelagic zone",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("Gulf", ignore_case = TRUE)) ~ "Gulf",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("sediment", ignore_case = TRUE)) ~ "Sediment",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("Sea bed", ignore_case = TRUE)) ~ "Benthic",
#     #str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("aquifer", ignore_case = TRUE)) ~ "Aquifer",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("AMALGI", ignore_case = TRUE)) ~ "Pelagic zone",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("P26_", ignore_case = TRUE)) ~ "Pelagic zone",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("P4_", ignore_case = TRUE)) ~ "Pelagic zone",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("JdFR", ignore_case = TRUE)) ~ "Bathypelagic",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("700m", ignore_case = TRUE)) ~ "Mesopelagic",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("1000", ignore_case = TRUE)) ~ "Bathypelagic",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("Malaspina", ignore_case = TRUE)) ~ "Bathypelagic",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("MSP", ignore_case = TRUE)) ~ "Bathypelagic",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("- MP", ignore_case = TRUE)) ~ "Bathypelagic",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("500_MG", ignore_case = TRUE)) ~ "Mesopelagic",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("200_MG", ignore_case = TRUE)) ~ "Mesopelagic",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("WHOI_OMZ", ignore_case = TRUE)) ~ "Pelagic zone",
#     str_detect(imgvr4_oceanic_metadata_filt$`PROJECT NAME`, regex("- LP-", ignore_case = TRUE)) ~ "Pelagic zone",
#     TRUE ~ `Ecosystem classification`  # leave everything else as-is
# ))
# #822,184 originally ... 574,750 w/o OMZ

# # Remove oceanic
# imgvr4_oceanic_metadata_filt_mod <- imgvr4_oceanic_metadata_filt_mod %>%
#   filter(`Ecosystem classification` != "Environmental;Aquatic;Marine;Oceanic")
# #546,664 after ... #299,230 w/o OMZ

# # See OMZs alone
# imgvr4_oceanic_metadata_OMZ <- imgvr4_oceanic_metadata_filt_mod %>%
#   filter(`Ecosystem classification` == "OMZ")
# #297,142 after ... #49,708 w/o OMZ

# imgvr4_oceanic_metadata_OMZ %>%
#   count(`PROJECT GOLD ID`, sort = TRUE)

# # Remove oceanic from OG so can recombine
# imgvr4_NoOceanic <- imgvr4_filt %>%
#   filter(`Ecosystem classification` != "Environmental;Aquatic;Marine;Oceanic")



# ###### Map PROJECT NAME to non oceanic to further filter if necessary ######
# # Make tax id col character rather than numeric
# imgvr4_NoOceanic$Taxon_oid <- as.character(imgvr4_NoOceanic$Taxon_oid)

# # Match taxid to table with GOLD ID
# imgvr4_NoOceanic <- gold_taxid %>%
#   dplyr::select(c("AP IMG TAXON ID", "AP PROJECT GOLD IDS")) %>%
#   right_join(imgvr4_NoOceanic, by = c("AP IMG TAXON ID" = "Taxon_oid"))

# # See how successful that was
# sum(is.na(imgvr4_NoOceanic$`AP PROJECT GOLD IDS`))
# #185 no match?

# # Match B and R cols in excel
# imgvr4_NoOceanic <- gold_metadata_filt %>%
#   right_join(imgvr4_NoOceanic, by = c("PROJECT GOLD ID" = "AP PROJECT GOLD IDS"))

# # See how successful that was
# sum(is.na(imgvr4_NoOceanic$`PROJECT NAME`))
# #10,549 with no project name



# # get smaller dfs for combining
# imgvr4_oceanic_metadata_filt_mod <- imgvr4_oceanic_metadata_filt_mod %>%
#   dplyr::select(c("UVIG", "Ecosystem classification"))

# imgvr4_NoOceanic <- imgvr4_NoOceanic %>%
#   dplyr::select(c("UVIG", "Ecosystem classification"))

# # Put the 2 together for plotting (non-oceanic and fixed oceanic)
# imgvr4_fixed <- bind_rows(imgvr4_oceanic_metadata_filt_mod, imgvr4_NoOceanic)


# ############################### Create plotting table ########################################

# # Group and sum for counts
# plot <- imgvr4_fixed %>%
#   group_by(`Ecosystem classification`) %>%
#   summarise(UVIG_count = n_distinct(UVIG))
# # Original sum here I have count 5 for Environmental;Aquatic;Marine; and count 2 for Environmental;Aquatic;Marine;Marine basin floor
# # Note that I am arbitrarily combining these with Environmental;Aquatic;Marine;Oceanic and Environmental;Aquatic;Marine;Benthic respectively

# # Add 5 and 2 to get rid of small categories that are indistinguishable
# plot <- plot %>%
#   mutate(`Ecosystem classification` = if_else(
#     `Ecosystem classification` == "Environmental;Aquatic;Marine;", # merging 5 marine viruses with oceanic
#     "Environmental;Aquatic;Marine;Pelagic",
#     `Ecosystem classification`)) %>%
#   group_by(`Ecosystem classification`) %>%
#   summarise(UVIG_count = sum(UVIG_count), .groups = "drop")

# plot <- plot %>%
#   mutate(`Ecosystem classification` = if_else(
#     `Ecosystem classification` == "Environmental;Aquatic;Marine;Marine basin floor",
#     "Environmental;Aquatic;Marine;Benthic", # merging 2 marine basin floor viruses with benthic
#     `Ecosystem classification`)) %>%
#   group_by(`Ecosystem classification`) %>%
#   summarise(UVIG_count = sum(UVIG_count), .groups = "drop")

# # Combine Pelagic and Pelagic zone
# plot <- plot %>%
#   mutate(`Ecosystem classification` = if_else(
#     `Ecosystem classification` == "Environmental;Aquatic;Marine;Pelagic",
#     "Environmental;Aquatic;Marine;Pelagic zone", # merging 471 pelagic viruses with pelagic zone
#     `Ecosystem classification`)) %>%
#   group_by(`Ecosystem classification`) %>%
#   summarise(UVIG_count = sum(UVIG_count), .groups = "drop")

# # Combine Hypoxic zone and OMZ
# plot <- plot %>%
#   mutate(`Ecosystem classification` = if_else(
#     `Ecosystem classification` == "Environmental;Aquatic;Marine;Hypoxic zone",
#     "Environmental;Aquatic;Marine;OMZ", # merging 9575 hypoxic zone viruses with OMZ
#     `Ecosystem classification`)) %>%
#   group_by(`Ecosystem classification`) %>%
#   summarise(UVIG_count = sum(UVIG_count), .groups = "drop")

# # Get the most specific name
# plot <- plot %>%
#   mutate(Specific_env = str_replace(`Ecosystem classification`, ";$", "")) %>%  # remove trailing ;
#   mutate(Specific_env = word(Specific_env, -1, sep = fixed(";"))) %>%           # take the last level
#   arrange(UVIG_count)

# # Get rid of redundant specific envs
# plot <- plot %>%
#   group_by(Specific_env) %>%
#   summarise(total_UVIGs = sum(UVIG_count, na.rm = TRUE))

# # how many viruses total in this?
# sum(plot$total_UVIGs)

# # add col of metadata for whether cultured or not
# plot <- plot %>%
#   mutate(Cultured = case_when(
#     Specific_env %in% c("Hydrothermal vents", "Ocean trench" , "Sediment", "Bathypelagic") ~ "yes",
#     TRUE ~ "no"
#   ))

################################ Lollipop plot #########################################

lolly <- ggplot(plot, aes(x = total_UVIGs, y = reorder(Specific_env, total_UVIGs))) +
  geom_segment(aes(x = 0, xend = total_UVIGs, yend = Specific_env), color = "#000000", linewidth = 0.7) +
  geom_point(size = 5, color = "#209bcf") +
  # geom_point(
  #   data = subset(plot, Cultured == "no"),
  #   aes(x = total_UVIGs + 0.05 * max(total_UVIGs)),
  #   shape = 17, size = 4, color = "red"
  # ) +
  labs(
    x = "UVIG Count",
    y = "Ecosystem",
    title = ""
  ) +
  theme_minimal(base_size = 15) +
  scale_x_continuous(expand = c(.01,0), labels = comma) +
  theme(
    panel.grid.major.y = element_blank(),  # remove horizontal gridlines
    panel.grid.minor.y = element_blank(),
    
    panel.grid.major.x = element_line(linetype = "dashed"),  # dashed vertical gridlines
    panel.grid.minor.x = element_blank()
  )
lolly

ggsave("/scratch/langwig/Projects/DeepSea_VirusReview/DeepSea_Vir/visuals/UVIGs_DeepSea_IMGVR.png", 
lolly, dpi = 500, bg = "transparent") #change height to 15 when mod'ing legend
