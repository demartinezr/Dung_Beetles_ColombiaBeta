setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata")

library(raster)
library(rgdal)
library(terra)
library(sf)
library(dplyr)

# extract America DEM
    divisions <- st_read("D:/Capas/America/countries/America_countries.shp")  # Shapefile con divisiones
    dem <- rast("D:/Capas/America/dem/srtm/srtm90m.tif")             # DEM
    divisions_vect <- vect(divisions)
    dem_crop <- crop(dem, divisions_vect)  
    dem_masked <- mask(dem_crop, divisions_vect)
#    writeRaster(dem_masked, "D:/Capas/World/america_srtm90m.tif", overwrite = TRUE)
#
# dataset
  records <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/GBIF_data/records_combined.rds")   # .rds file with GBIF records
  dem <- rast("D:/Capas/America/dem/srtm/america_srtm90m.tif")           # DEM file (tif)
#
# get elevation from strm 90m DEM
  records <- records[complete.cases(records$decimalLongitude, records$decimalLatitude), ]
  records_sf <- st_as_sf(records, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  elevation <- extract(dem, vect(records_sf))
  records$elevation <- elevation$elevation
  saveRDS(records, file="./GBIF_data/records_combined_ele.rds")
  registros <- readRDS("./GBIF_data/records_combined_ele.rds")
  # <- as.data.frame(registros@data)
  sp_ele <- registros %>%
    group_by(scientificName1) %>%
    summarise(
      LimiteInferior = min(elevation, na.rm = TRUE),
      LimiteSuperior = max(elevation, na.rm = TRUE),
      Promedio = mean(elevation, na.rm = TRUE),
      ConteoRegistros = n()
    ) %>%
    mutate(LimiteInferior = ifelse(LimiteInferior < 0, 0, LimiteInferior)) %>%
    ungroup()
  
  write.csv(sp_ele, "sp_ele_strm90m.csv", row.names = FALSE)
  
