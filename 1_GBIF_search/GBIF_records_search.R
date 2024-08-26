# Directory
    setwd("D:/Doctorado/Tesis/GBIF")
# Packages
    library(readxl)
    library(rgbif)
    library(dplyr)
    library(ggplot2)
    library(rnaturalearth)
    library(rnaturalearthdata)
    library(sf)

# Dung beetles Dataset #
    IBDdata <- read_excel("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/Scarabaeinae_database_2024.xlsx", sheet="Scarabaeinae_database_2024")
#
# List of species and subspecies #
    species <- data.frame(subset(IBDdata, taxonRank == "Species"))
    species <- sort(unique(species$scientificName)) 
    species <- gsub("_", " ", species)
    subspecies <- data.frame(subset(IBDdata, taxonRank == "Subspecies"))
    subspecies <- sort(unique(subspecies$scientificName)) 
    subspecies <- gsub("_", " ", subspecies)
    species <- c(species, subspecies)
#
# GBIF search #
    result <- occ_data(scientificName = species, hasCoordinate = TRUE, synonym = TRUE)
    data_list <- lapply(result, function(x) x$data)
    records <- bind_rows(data_list, .id = "scientificName1")
    #names(records)
    #write.table(records, "GBIF_2024-02-13.txt", row.names = FALSE, sep=";")
#
# Check records #
#
# subset records of Deltochilum carinatum
    d_car <- data.frame(subset(records, scientificName1 == "Deltochilum carinatum"))
    View(d_car)
#
# base maps #
    south_america <- ne_countries(scale = "medium", continent = "South America", returnclass = "sf")
    central_america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
    central_america <- central_america[central_america$name %in% c("Belize", "Costa Rica", "El Salvador", "Mexico", "Guatemala", "Honduras", "Nicaragua", "Panama"), ]
    ecoreg <- st_read(dsn = "D:/Capas/World/wwf/terr-ecoregions-TNC", layer = "tnc_terr_ecoregions") #Pais
    ecoreg <- st_transform(ecoreg, CRS("+proj=longlat +datum=WGS84"))
    col <- ne_states(country = "colombia", returnclass = "sf")
# plot 
    ggplot() +
     geom_sf(data = south_america, fill = "lightgrey", color = "black") +
     geom_sf(data = central_america, fill = "lightgrey", color = "black") +
     geom_point(data = d_car, aes(x = decimalLongitude, y = decimalLatitude), color = "blue")
