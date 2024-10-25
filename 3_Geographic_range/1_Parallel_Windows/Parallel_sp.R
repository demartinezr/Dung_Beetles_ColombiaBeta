# Function for obtaining geographic ranges based on GBIF records, ecoregions, and 
# altitudinal ranges (alos palsar DEM 30m) for species.
#
  setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata")
#
# R packages
  library(raster)
  library(rgdal)
  library(sf)
  library(rgeos)
  library(terra)
  library(dplyr)
  library(doParallel)
  library(foreach)
  
  # Load the data and layers
  records <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/GBIF_data/records_combined_ele.rds")
  dem <- rast("D:/Capas/America/dem/srtm/col_srtm90m.tif")
  divisions <- st_read("D:/Capas/Colombia/Colombia/COL_adm0.shp")
  # Get unique species and limit to 243
  species_subset <- sort(unique(records$scientificName1))[1:5]
  
  # Define the function to process each species
  process_species <- function(species, records, divisions, dem) {
    cat("Processing species:", species, "\n")
    
    each_species <- records[records$scientificName1 == species, ]
    coordinates <- each_species[, c("decimalLongitude", "decimalLatitude")]
    coordinates <- na.omit(coordinates)  # Remove rows with NA values
    coordinates <- coordinates[!duplicated(coordinates), ]  # Remove duplicate coordinates
    
    if (nrow(coordinates) == 0) {
      cat("No records for:", species, "\n")
      return(NULL)  # Return NULL if there are no records
    }
    
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
    
    polygon_sf <- st_transform(polygon_sf, crs = 4326)
    polygon_col <- st_intersection(polygon_sf, divisions)
    
    if (nrow(polygon_col) > 0) {
      dem_contour <- crop(dem, polygon_col)
      dem_contour <- mask(dem_contour, polygon_col)
      
      max_elevation <- max(each_species$elevation, na.rm = TRUE) + 100
      min_elevation <- max(min(each_species$elevation, na.rm = TRUE) - 100, 0)
      
      adjusted_altitude_raster <- dem_contour >= min_elevation & dem_contour <= max_elevation
      contours <- terra::as.polygons(adjusted_altitude_raster, dissolve = TRUE)
      contours <- contours[contours$elevation == 1, ]
      
      file_name <- paste0(gsub(" ", "_", species), ".shp")
      file_path <- file.path("D:/Doctorado/Tesis/GBIF/contorno", file_name)
      writeVector(contours, file_path, filetype = "ESRI Shapefile", overwrite = TRUE)
      cat("Contour saved as:", file_path, "\n")
    } else {
      cat("No overlapping polygons found for:", species, "\n")
    }
  # Clean up memory after processing
  rm(each_species, coordinates, points, polygon_sf, polygon_col, dem_contour, adjusted_altitude_raster, contours)
  gc()  # Call garbage collector to free memory
  }
  # Set up parallel backend to use 4 cores
  cl <- makeCluster(4)
  
  # Export required variables to the cluster
  clusterExport(cl, varlist = c("records", "divisions", "dem", "st_as_sf", "st_buffer", "st_transform", "st_union", "st_convex_hull", "crop", "mask", "terra::as.polygons", "writeVector"))
  
  registerDoParallel(cl)
  
  # Parallel processing
  species_results <- foreach(species = species_subset, .combine = 'c', 
                             .packages = c('sf', 'terra', 'dplyr')) %dopar% {
                               process_species(species, records, divisions, dem)
                               gc()  # Force garbage collection to free memory after each iteration
                               return(result)
                               }
  
  stopCluster(cl)
  