# ------------------------------------------------------------------------------
# Elevation range and standardized elevation for dung beetle species
#
# This script compiles elevation records for dung beetle species and morphospecies
# from multiple sources, including GBIF occurrences, SIB-Colombia and entomological
# collections, project sampling data, and additional literature-based records.
#
# Elevation values were extracted from the ALOS PALSAR Digital Elevation Model
# (30 m resolution) and used to estimate species-level elevation ranges.
# For each taxon, minimum, maximum, and mean elevation were calculated from
# unique geographic records. When only a single record was available, the
# elevation range was expanded by ±100 m to account for potential measurement
# uncertainty and limited sampling.
#
# The resulting elevation ranges were then standardized to produce a scaled
# elevation variable (−1 = minimum elevation, 0 = mean elevation, 1 = maximum
# elevation). This standardized elevation gradient is used as an environmental
# predictor in the Bayesian hierarchical model of dung beetle abundance,
# allowing the estimation of species-specific responses to elevation and
# land-use change.
# ------------------------------------------------------------------------------
  setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata")
#
# R Packages
  library(readxl)
  library(dplyr)
#
# Load species list with elevation from alos palsar dem 30m 
  # species records from GBIF
    sp <- readRDS("./elevation_range/registros_ele.rds")
    sp <- sp %>%
      group_by(scientificName1) %>%
      distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
      ungroup()
  #
  # morphospecies verified with IAvH code "H" from SIB-Colombia and collections
    msp <- readRDS("./elevation_range/morphoIAvH_SIB_ele.rds")
    msp <- mutate(msp, identificationRemarks = ifelse(identificationRemarks=="Canthidium_sp_13H", "Canthidium sp. 13H", identificationRemarks))
    msp <- mutate(msp, identificationRemarks = ifelse(identificationRemarks=="Canthon_sp_11H", "Canthon sp. 11H", identificationRemarks))
    msp <- mutate(msp, identificationRemarks = ifelse(identificationRemarks=="Canthon_sp_18H", "Canthon sp. 18H", identificationRemarks))
    msp <- mutate(msp, identificationRemarks = ifelse(identificationRemarks=="Canthon_sp_19H", "Canthon sp. 19H", identificationRemarks))
    msp <- mutate(msp, identificationRemarks = ifelse(identificationRemarks=="Canthon_sp_22H", "Canthon sp. 22H", identificationRemarks))
    msp <- mutate(msp, identificationRemarks = ifelse(identificationRemarks=="Deltochilum_sp_18H", "Deltochilum sp. 18H", identificationRemarks))
    msp <- mutate(msp, identificationRemarks = ifelse(identificationRemarks=="Onthophagus_sp_08H", "Onthophagus sp. 08H", identificationRemarks))
    msp <- mutate(msp, identificationRemarks = ifelse(identificationRemarks=="Uroxys_sp_03H", "Uroxys sp. 03H", identificationRemarks))
    msp <- mutate(msp, identificationRemarks = ifelse(identificationRemarks=="Uroxys_sp_06H", "Uroxys sp. 06H", identificationRemarks))
    msp <- mutate(msp, identificationRemarks = ifelse(identificationRemarks=="Uroxys_sp_08H", "Uroxys sp. 08H", identificationRemarks))
    msp <- msp %>%
      group_by(identificationRemarks) %>%
      distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
      ungroup()
  #
  # morphospecies from project, using our data
    uvmsp <- read_excel("./abundance/Scarabaeinae_database_2024.xlsx", sheet="Scarabaeinae_database_2024")
    uvmsp <- data.frame(subset(uvmsp, taxonRank == "Unverified morphospecies")) 
    uvmsp <- uvmsp %>%
    rename(decimalLatitude = lat_all_points,
          decimalLongitude = lon_all_points)
    uvmsp$decimalLatitude <- as.numeric(uvmsp$decimalLatitude)
    uvmsp$decimalLongitude <- as.numeric(uvmsp$decimalLongitude)
    uvmsp$elev_ALOS_all_points <- as.numeric(uvmsp$elev_ALOS_all_points)
    uvmsp <- uvmsp %>%
      group_by(scientificName) %>%
      distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
      ungroup()
  #
  # Other species/morphospecies with not records in GBIF and SIB Colombia
    sp_proj <- read_excel("./abundance/Scarabaeinae_database_2024.xlsx", sheet="Scarabaeinae_database_2024")
    sp_proj <- data.frame(subset(sp_proj, scientificName %in% c("Gromphas_aeruginosa", 
      "Ontherus_sulcator", "Uroxys_pygmaeus", "Ateuchus_cracicus", "Canthidium_sp._20H", 
      "Coprophanaeus_edmondsi", "Cryptocanthon_campbellorum", "Cryptocanthon_mailinae", 
      "Deltochilum_sp_.21H", "Deltochilum_sp._22H", "Deltochilum_sp._27H", "Deltochilum_sp._28H", 
      "Deltochilum_sp_.29H", "Deltochilum_sp._31H", "Deltochilum_sp._32H", "Deltochilum_sp._33H", 
      "Deltochilum_sp_.35H", "Deltochilum_sp._38H", "Deltochilum_sp._39H", "Deltochilum_sp._48H", 
      "Deltochilum_sp_.54H", "Deltochilum_sp._55H", "Ontherus_sanctamartae")))
    sp_proj <- sp_proj %>%
      rename(decimalLatitude = lat_all_points,
            decimalLongitude = lon_all_points)
    sp_proj$decimalLatitude <- as.numeric(sp_proj$decimalLatitude)
    sp_proj$decimalLongitude <- as.numeric(sp_proj$decimalLongitude)
    sp_proj$elev_ALOS_all_points <- as.numeric(sp_proj$elev_ALOS_all_points)
    sp_proj <- sp_proj %>%
      group_by(scientificName) %>%
      distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
      ungroup()
  # species with no records from literature, upload D. compressicollis
    spp <- readRDS("./elevation_range/varios_ele.rds")
    spp <- spp %>%
      group_by(scientificName) %>%
      distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE) %>%
      ungroup()
  # select columns for datasets
    sp <- sp[c("scientificName1", "alos_ele")]
    sp <- sp %>% rename(scientificName = scientificName1)
    msp <- msp[c("identificationRemarks", "alos_ele")]
    msp <- msp %>% rename(scientificName = identificationRemarks)
    uvmsp <- uvmsp[c("scientificName", "elev_ALOS_all_points")]
    uvmsp <- uvmsp %>% rename(alos_ele = elev_ALOS_all_points)
    sp_proj <- sp_proj[c("scientificName", "elev_ALOS_all_points")]
    sp_proj <- sp_proj %>% rename(alos_ele = elev_ALOS_all_points)
    spp <- spp[c("scientificName", "alos_ele")]
  #
    # Combine all elevation records into a unified dataset
    species_join <- bind_rows(sp, msp, uvmsp, sp_proj, spp)
#
    # Estimate species-level elevation range
    # If only one record is available, extend the range by ±100 m
  elevation_range <- species_join %>%
    group_by(scientificName) %>%
    summarise(
      count = n(),
      min_elevation = ifelse(n() == 1, pmax(alos_ele - 100, 0), min(alos_ele, na.rm = TRUE)),
      max_elevation = ifelse(n() == 1, alos_ele + 100, max(alos_ele, na.rm = TRUE)),
      avg_elevation = ifelse(n() == 1, alos_ele, mean(alos_ele, na.rm = TRUE))
    ) %>%
    ungroup()
# write.table(elevation_range, file="./elevation_range/elevation_range_2024.txt", sep=",")
#
  # Standardize elevation to create a relative elevation gradient
  # (-1 = minimum elevation, 0 = mean elevation, 1 = maximum elevation)
  # This standardized variable is used as the elevation predictor in the Bayesian model
    elevation_standardized <- elevation_range %>%
    mutate(
      min_elevation = -1,
      max_elevation = 1,
      avg_elevation = 0
    )
    