setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata")

library(raster)
library(rgdal)
library(sf)
library(dplyr)

# dataset
  records <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/GBIF_data/records_combined.rds")   # .rds file with GBIF records
  dem <- rast("C:/Users/Dell-PC/Dropbox/CO_DBdata/SIG/coldem30/alos_dem_col.tif")           # DEM file (tif)
#
# get elevation from Alos palsar 30m DEM
  records <- records[complete.cases(records$decimalLongitude, records$decimalLatitude), ]
  records_sf <- st_as_sf(records, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  elevation <- extract(dem, vect(records_sf))
  records$alos_ele <- elevation$alos_dem_col
  saveRDS(records, file="./GBIF_data/records_combined_ele.rds")
  registros <- readRDS("./GBIF_data/records_combined_ele.rds")
  # <- as.data.frame(registros@data)
  sp_ele <- registros %>%
    group_by(scientificName1) %>%
    summarise(
      LimiteInferior = min(alos_ele, na.rm = TRUE),
      LimiteSuperior = max(alos_ele, na.rm = TRUE),
      Promedio = mean(alos_ele, na.rm = TRUE)
    ) %>%
  ungroup() 
