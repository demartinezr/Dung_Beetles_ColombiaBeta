# Protocol to obtain geographic ranges from GBIF records, convex hulls, and 
# altitudinal ranges (DEM)
#
# Libraries
  library(raster)
  library(dplyr)
  library(ggplot2)
  library(rgdal)
  library(sf)
  library(rgeos)
  library(terra)
#
  # dataset and Colombia shapefile and DEM STRM 90m
  records <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/GBIF_data/records_combined_ele.rds")
  dem <- rast("D:/Capas/America/dem/srtm/col_srtm90m.tif")
  divisions <- st_read("D:/Capas/Colombia/Colombia/COL_adm0.shp")
#
# Geographic ranges for a single species
#
# Select the species
  species <- sort(unique(records$scientificName1))[52]
  cat("Processing species:", species, "\n")
  each_species <- records[records$scientificName1 == species, ]
  # Get coordinates for the selected species
  coordinates <- each_species[, c("decimalLongitude", "decimalLatitude")]
  coordinates <- na.omit(coordinates)  # Remove rows with NA values
  coordinates <- coordinates[!duplicated(coordinates), ]  # Remove duplicate coordinates
  # Check if there are records for the species
  if (nrow(coordinates) == 0) {
    cat("No records for:", species, "\n")
  } else {
    # Convert coordinates to a spatial object
    points <- st_as_sf(coordinates, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    # Create a polygon based on the number of records
    if (nrow(coordinates) == 1) {
      polygon_sf <- st_buffer(st_transform(points, crs = 3857), dist = 10000)
    } else if (nrow(coordinates) == 2) {
      line <- st_union(points) %>% st_cast("LINESTRING")
      polygon_sf <- st_buffer(st_transform(line, crs = 3857), dist = 10000)
    } else {
      polygon_sf <- st_convex_hull(st_union(points))
      polygon_sf <- st_buffer(st_transform(polygon_sf, crs = 3857), dist = 10000) 
    }
  # Adjust the polygon to the divisions of Colombia
    polygon_sf <- st_transform(polygon_sf, crs = 4326)
    polygon_col <- st_intersection(polygon_sf, divisions)
    polygon_col <- vect(polygon_col)
    if (nrow(polygon_col) > 0) {
      # Crop and mask the DEM with the adjusted polygon
      dem_contour <- crop(dem, polygon_col)
      dem_contour <- mask(dem_contour, polygon_col)
      # Define the altitude limits
      max_elevation <- max(each_species$elevation, na.rm = TRUE) + 100
      min_elevation <- max(min(each_species$elevation, na.rm = TRUE) - 100, 0)
      # Create a logical raster with cells within the altitude range
      adjusted_altitude_raster <- dem_contour >= min_elevation & dem_contour <= max_elevation
      # Get the contour of the adjusted raster
      contours <- terra::as.polygons(adjusted_altitude_raster, dissolve = TRUE)
      contours <- contours[contours$elevation == 1, ]
      # Save the resulting contour
      file_name <- paste0(gsub(" ", "_", species), ".shp")
      file_path <- file.path("D:/Doctorado/Tesis/GBIF/contorno", file_name)
      writeVector(contours, file_path, filetype = "ESRI Shapefile", overwrite=TRUE)
      
      cat("Contour saved as:", file_path, "\n")
    } else {
      cat("No overlapping polygons found for:", species, "\n")
    }
  }
  
###############################################################################
# Get geographic ranges for a multiple species
#  
  library(raster)
  library(dplyr)
  library(ggplot2)
  library(rgdal)
  library(sf)
  library(rgeos)
  library(terra)
  #
# Dataset, Colombia shapefile and DEM SRTM 90m
  records <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/GBIF_data/records_combined_ele.rds")
  dem <- rast("D:/Capas/America/dem/srtm/col_srtm90m.tif")
  divisions <- st_read("D:/Capas/Colombia/Colombia/COL_adm0.shp")
# function for multiples species
  #
  process_species <- function(species) {
    cat("Processing species:", species, "\n")
    each_species <- records[records$scientificName1 == species, ]
    # Get coordinates pof species i
    coordinates <- each_species[, c("decimalLongitude", "decimalLatitude")]
    coordinates <- na.omit(coordinates)
    coordinates <- coordinates[!duplicated(coordinates), ]
    # Check if there are records for the species
    if (nrow(coordinates) == 0) {
      cat("No records for:", species, "\n")
      return(NULL)  # Return NULL if there are no records
    }
    # Convert coordinates to a spatial object
    points <- st_as_sf(coordinates, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    # Create a polygon based on the number of records
    if (nrow(coordinates) == 1) {
      polygon_sf <- st_buffer(st_transform(points, crs = 3857), dist = 50000)
    } else if (nrow(coordinates) == 2) {
      line <- st_union(points) %>% st_cast("LINESTRING")
      polygon_sf <- st_buffer(st_transform(line, crs = 3857), dist = 50000)
    } else {
      polygon_sf <- st_convex_hull(st_union(points))
      polygon_sf <- st_buffer(st_transform(polygon_sf, crs = 3857), dist = 50000) 
    }
    # Adjust the polygon to the divisions of Colombia
    polygon_sf <- st_transform(polygon_sf, crs = 4326)
    polygon_col <- st_intersection(polygon_sf, divisions)
    polygon_col <- vect(polygon_col)
    if (nrow(polygon_col) > 0) {
      # Crop and mask the DEM with the adjusted polygon
      dem_contour <- crop(dem, polygon_col)
      dem_contour <- mask(dem_contour, polygon_col)
      # Define the altitude limits
      max_elevation <- max(each_species$elevation, na.rm = TRUE) + 100
      min_elevation <- max(min(each_species$elevation, na.rm = TRUE) - 100, 0)
      # Create a logical raster with cells within the altitude range
      adjusted_altitude_raster <- dem_contour >= min_elevation & dem_contour <= max_elevation
      # Get the contour of the adjusted raster
      contours <- terra::as.polygons(adjusted_altitude_raster, dissolve = TRUE)
      contours <- contours[contours$elevation == 1, ]
      if (nrow(contours) > 0) {
        contours$scientificName <- gsub(" ", "_", species)
      # Save the resulting contour
      file_name <- paste0(gsub(" ", "_", species), ".shp")
      file_path <- file.path("D:/Doctorado/Tesis/GBIF/contorno", file_name)
      writeVector(contours, file_path, filetype = "ESRI Shapefile", overwrite = TRUE)
      cat("Contour saved as:", file_path, "\n")
    } else {
      cat("No overlapping polygons found for:", species, "\n")
    }
    }
   # Release memory
    rm(each_species, coordinates, points, polygon_sf, polygon_col, dem_contour, adjusted_altitude_raster, contours)
    gc()  # Call garbage collector
  }
  # Process each species in the species_subset
  for (species in species_subset) {
    process_species(species)
  }
  
# combined the 243 shapefiles
  shapefiles_list <- list.files("D:/Doctorado/Tesis/GBIF/contorno", pattern = "*.shp", full.names = TRUE)
  shapefiles_combined <- lapply(shapefiles_list, st_read) %>% do.call(rbind, .)
  saveRDS(shapefiles_combined, "C:/Users/Dell-PC/Dropbox/CO_DBdata/geographic_range/geographic_range.rds")
  #