setwd("D:/Doctorado/Tesis/GBIF")

library(raster)
#library(rgbif)
library(readxl)
library(dplyr)
library(maps)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(rgdal)
library(sp)
library(sf)
library(rgeos)

#################################### Data #####################################
IBDdata <- read_excel("C:/Users/Dell-PC/Dropbox/CO_DBdata/Scarabaeinae_database_2024.xlsx", sheet="Scarabaeinae_database")
####################### List of species and subspecies ########################
species <- data.frame(subset(IBDdata, taxonRank == "species"))
species <- sort(unique(species$scientificName)) 
species <- gsub("_", " ", species)
subspecies <- data.frame(subset(IBDdata, taxonRank == "Subspecies"))
subspecies <- sort(unique(subspecies$scientificName)) 
subspecies <- gsub("_", " ", subspecies)
species <- c(species, subspecies)
############################## GBIF search ####################################
resultado <- occ_data(scientificName = species, hasCoordinate = TRUE, synonym = TRUE)
#names(registros)
#names(registros[["Anisocanthon villosus"[1]]]$data)
datos_lista <- lapply(resultado, function(x) x$data)
registros <- bind_rows(datos_lista, .id = "scientificName1")
#registros <- do.call(cbind, registros) #exportar a txt
#registros <- as.data.frame(registros) #exportar a txt
#write.table(registros, "GBIF_2024-02-13.txt", row.names = FALSE, sep=";")
############################ plots ####################################
### subset records of Deltochilum carinatum
d_car <- data.frame(subset(registros, scientificName1 == "Deltochilum carinatum"))
View(d_car)
#################################### base maps #################################
south_america <- ne_countries(scale = "medium", continent = "South America", returnclass = "sf")
central_america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
central_america <- central_america[central_america$name %in% c("Belize", "Costa Rica", "El Salvador", "Mexico", "Guatemala", "Honduras", "Nicaragua", "Panama"), ]
ecoreg <- readOGR(dsn = "D:/Users/DELL/Documents/Capas/World/wwf/terr-ecoregions-TNC", layer = "tnc_terr_ecoregions") #Pais
ecoreg <- spTransform(ecoreg, CRS(proj4string(poligono_sp)))
col <- ne_states(country = "colombia", returnclass = "sf")
#shp <- raster::getData('GADM', country = 'COL', level = 1)
# plot 
ggplot() +
  geom_sf(data = south_america, fill = "lightgrey", color = "black") +
  geom_sf(data = central_america, fill = "lightgrey", color = "black") +
#  geom_sf(data = sol_com, fill = "green", color = "black") +
  geom_point(data = d_car, aes(x = decimalLongitude, y = decimalLatitude), color = "blue")
