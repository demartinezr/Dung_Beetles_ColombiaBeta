# This script simulates the posterior occupancy probability for each species-point across Colombia at 2 km resolution
# for one posterior iteration

library(sf)
library(reticulate)
library(raster)

`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

##### Get elevation raster across Colombia #####
# Set up GEE session
use_condaenv('rgee', conda = "auto", required = TRUE) # point reticulate to the conda environment created in GEE_setup.sh
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()            # Trigger the authentication

np <- import("numpy")       # Import Numpy        needed for converting gee raster to R raster object
pd <- import("pandas")      # Import Pandas       ditto the above

# Get elevations for Colombia
countries <- ee$FeatureCollection('USDOS/LSIB_SIMPLE/2017')
roi <- countries$filterMetadata('country_na', 'equals', 'Colombia')  
DEM_full <- ee$Image("JAXA/ALOS/AW3D30/V2_2")
DEM <- DEM_full$select(list("AVE_DSM"),list("elevation"))
latlng <- ee$Image$pixelLonLat()$addBands(DEM)
latlng <- latlng$reduceRegion(reducer = ee$Reducer$toList(),
                              geometry = roi,
                              maxPixels = 10^9,
                              scale=2000)
# Convert to arrays
lats <- np$array((ee$Array(latlng$get("latitude"))$getInfo()))
lngs <- np$array((ee$Array(latlng$get("longitude"))$getInfo()))
elevs <- np$array((ee$Array(latlng$get("elevation"))$getInfo()))
# Convert to elevation raster
elevation <- data.frame(x = lngs, y = lats, elevation = elevs)
raster_elev <- raster::rasterFromXYZ(elevation,crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# Reproject
raster_elev_AEA <- raster::projectRaster(raster_elev, crs = sp::CRS(AEAstring))
# Save raster
# raster::writeRaster(raster_elev, "D:/Capas/America/dem/elev_raster/raster_elev.grd", overwrite = T)
# raster::writeRaster(raster_elev_AEA, "D:/Capas/America/dem/elev_raster/raster_elev_AEA.grd", overwrite = T)

raster_elev <- raster::raster("D:/Capas/America/dem/elev_raster/raster_elev.grd")
raster_elev_AEA <- raster::raster("D:/Capas/America/dem/elev_raster/raster_elev_AEA.grd")

##### Rasterize the goegraphic range maps
#geo_range_shape <- list.files("D:/Doctorado/Tesis/GBIF/contorno", pattern = "\\.shp$", full.names = TRUE)
#geo_range_list <- lapply(geo_range_shape, st_read)
#selected_shapefiles <- do.call(rbind, geo_range_list)
#selected_shapefiles <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/geographic_range/geographic_range.rds")
#saveRDS(geo_range_list, "D:/Capas/America/geo_ranges/geo_range_shapes/geo_range_shapes.rds")
geo_range_shapes <- readRDS("D:/Capas/America/geo_ranges/geo_range_shapes/geo_range_shapes.rds")
geo_range_shapes_AEA <- lapply(geo_range_shapes, function(shp) {
  st_transform(shp, crs(raster_elev_AEA))})

pb <- txtProgressBar(min = 0, max = length(geo_range_shapes_AEA), initial = 0, style = 3) 
for(i in 1:length(geo_range_shapes_AEA)){
  setTxtProgressBar(pb,i)
  geo_range_raster <- fasterize::fasterize(st_sf(st_union(geo_range_shapes_AEA[[i]])), raster_elev_AEA)
  raster::writeRaster(geo_range_raster, paste0('D:/Capas/America/geo_ranges/geo_range_rasters/geo_range_rasters', species_list[i], '_buffered.grd'), overwrite = T)
}

#### Distance-to-range for each species ####
library(brms)
db_mod_abundance <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod_abundance.rds")
db7 <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db7_combined.rds")
raster_elev_AEA <- raster::raster("D:/Capas/America/dem/elev_raster/raster_elev_AEA.grd")
species_list <- unique(db7$scientificName)

dtr_mean <- mean(db7$distance_from_range)
dtr_sd <- sd(db7$distance_from_range)

points_for_dist <- raster_elev_AEA
mask_raster <- raster("D:/Capas/America/dem/mask/mask.grd")
points_for_dist[mask_raster == 1] <- NA

pfd <- raster::as.data.frame(points_for_dist, xy=T)
pfd <- pfd[!is.na(pfd$elevation),]

pfd2 <- coordinates(raster_elev_AEA, spatial = T)

all.equal(pfd2@coords[,1], pfd$x)
all.equal(pfd2@coords[,2], pfd$y)

points <- st_as_sf(pfd2[!is.na(pfd$elevation), ]) 
#points <- st_as_sf(pfd, coords = c("x", "y"), crs = st_crs(raster_elev_AEA))

coords <- st_coordinates(points)

sp_list <- species_list
names(geo_range_shapes) <- species_list
geo_range_shapes_AEA <- lapply(geo_range_shapes, function(shp) {
  st_transform(shp, crs(raster_elev_AEA))})

pb <- txtProgressBar(min = 0, max = length(sp_list), initial = 0, style = 3) 
for(i in 1:length(sp_list)){
  setTxtProgressBar(pb,i)
  sp <- sp_list[i]
  geo_range_polygon <- st_union(geo_range_shapes_AEA[[sp]])
  distances <- as.numeric(st_distance(points, st_cast( geo_range_polygon, to = "MULTILINESTRING")))*(2*(as.numeric(st_distance(points, geo_range_polygon))>0) - 1) # The second part gives positive distances for outside-of-range and negative distances for in-range.  Turns out that as.numeric(st_distance(points, ayerbe))>0) is much faster than !st_within(points, ayerbe) 
  transformed_distances <- boot::inv.logit(4.7*(distances - dtr_mean)/dtr_sd)
  td_df <- cbind(coords, transformed_distances)
  range_dist_raster_i <- raster::rasterFromXYZ(td_df,crs=AEAstring)
  raster::writeRaster(range_dist_raster_i, paste0('D:/Capas/America/geo_ranges/transformed_distance/', sp_list[i], '.grd'), overwrite = F)
}

## predict abundance across Colombia ##

raster_elev_AEA <- raster::raster("D:/Capas/America/dem/elev_raster/raster_elev_AEA.grd")
elev_df <- raster::as.data.frame(raster_elev_AEA, xy = T) 
elev_df$cell_id <- 1:nrow(elev_df)
# names(elev_df)[3] <- "elevation"
lambda <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/Analysis/get_posterior/posterior_lambda.rds")
db_mod_abundance <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod_abundance.rds")
draws <- posterior::as_draws_df(db_mod_abundance)
# draws <- draws[1:2000, ]
db7 <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db7_combined.rds")




