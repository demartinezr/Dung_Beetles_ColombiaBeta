# This script generates geographic range polygons for dung beetle species using 
# GBIF and biological collection records, regional geographic divisions, and elevation 
# limits derived from a 30-m resolution DEM. For each species, the workflow constructs 
# a geometry based on the number of available occurrence points: a 10-km buffer around 
# a single point, a 10-km buffered line when only two records are present, and a 10-km 
# buffered convex hull when three or more records are available. These preliminary 
# ranges are clipped to the corresponding geographic divisions and further refined by 
# applying species-specific altitudinal limits, defined as the minimum and maximum 
# recorded elevations expanded by ±100 meters to account for local uncertainty. 
# The resulting elevation-filtered polygons represent the estimated geographic range 
# and are exported as ESRI Shapefiles.

  setwd("C:/Users/PC/Dropbox/CO_DBdata")
#
# R packages
  library(raster)
  library(sf)
  library(terra)
  library(dplyr)
  library(doParallel)
  library(foreach)
  
  # Load the data and layers
  records <- readRDS("./GBIF_data/records_combined_ele.rds")
  dem <- rast("F:/Capas/America/dem/elev_raster/raster_elev.grd")
  source("F:/repositorio/Dung_Beetles_ColombiaBeta/3_Geographic_range/mainland/mainland.R")
  
  records_summary <- records %>%
    count(scientificName1, source) %>%
    pivot_wider(names_from=source,
                values_from = n,
                values_fill = 0) %>%
    mutate(total = rowSums(across(where(is.numeric))))
  
  write.table(records_summary, file="records_summary.txt",
              sep=",", row.names = FALSE, quote=FALSE)
  
  # Get unique species and limit to 243
  species_subset <- sort(unique(records$scientificName1))[1:2]
  
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
  