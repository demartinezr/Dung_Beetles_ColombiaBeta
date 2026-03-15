###############################################################################
# Geographic range estimation from occurrence records
#
# This script estimates geographic ranges for species using occurrence records
# obtained in 2_GBIF_clean module.
# Following the approach described in the main manuscript and
# Supplemental Methods, species ranges are approximated using the extent of
# occurrence (EOO) derived from spatial records.
#
# Occurrence coordinates are first filtered to retain unique spatial records.
# Depending on the number of available records per species, geographic ranges
# are estimated using different spatial envelopes: a buffered point for single
# records, a buffered line connecting two records, or a convex hull constructed
# from three or more records. All geometries are buffered to account for
# georeferencing uncertainty and potential unsampled areas near known records.
#
# The resulting polygons are intersected with the political boundary of
# Colombia to restrict species ranges to the study area. Elevational limits are
# then applied using a digital elevation model (DEM), where the observed
# minimum and maximum elevations from occurrence records are expanded by
# ±100 m to accommodate minor measurement errors and small range extensions.
# The DEM is masked using these elevation limits and converted to polygons,
# representing the potential geographic range of each species within Colombia.
#
# Individual species ranges are exported as shapefiles and subsequently merged
# into a single dataset for downstream analyses.
###############################################################################
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
  dem <- rast("D:/Capas/America/dem/ALOS/col_alos30m.tif")
  divisions <- st_read("D:/Capas/Colombia/Colombia/COL_adm0.shp")
# Function to estimate the geographic range for a given species
  #
  process_species <- function(species) {
    cat("Processing species:", species, "\n")
    each_species <- records[records$scientificName1 == species, ]
    # Extract occurrence records for the focal species
    coordinates <- each_species[, c("decimalLongitude", "decimalLatitude")]
    coordinates <- na.omit(coordinates)
    coordinates <- coordinates[!duplicated(coordinates), ]
    # Check if there are records for the species
    if (nrow(coordinates) == 0) {
      cat("No records for:", species, "\n")
      return(NULL)  # Return NULL if there are no records
    }
    # Convert occurrence coordinates to spatial points
    points <- st_as_sf(coordinates, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    # Generate a preliminary geographic envelope depending on the number of records
    # Single record: buffered point representing a minimal spatial extent
    if (nrow(coordinates) == 1) {
      polygon_sf <- st_buffer(st_transform(points, crs = 3857), dist = 50000)
      # Two records: buffered line connecting both occurrences
    } else if (nrow(coordinates) == 2) {
      line <- st_union(points) %>% st_cast("LINESTRING")
      polygon_sf <- st_buffer(st_transform(line, crs = 3857), dist = 50000)
      # Three or more records: convex hull representing the extent of occurrence (EOO)
    } else {
      polygon_sf <- st_convex_hull(st_union(points))
      polygon_sf <- st_buffer(st_transform(polygon_sf, crs = 3857), dist = 50000) 
    }
    # Restrict the preliminary polygon to the Colombian national boundary
    polygon_sf <- st_transform(polygon_sf, crs = 4326)
    polygon_col <- st_intersection(polygon_sf, divisions)
    polygon_col <- vect(polygon_col)
    if (nrow(polygon_col) > 0) {
      # Crop and mask the DEM using the species geographic envelope
      dem_contour <- crop(dem, polygon_col)
      dem_contour <- mask(dem_contour, polygon_col)
      # Identify DEM cells falling within the adjusted elevation range
      max_elevation <- max(each_species$elevation, na.rm = TRUE) + 100
      min_elevation <- max(min(each_species$elevation, na.rm = TRUE) - 100, 0)
      # Convert suitable elevation cells to polygons representing potential habitat
      adjusted_altitude_raster <- dem_contour >= min_elevation & dem_contour <= max_elevation
      # Get the contour of the adjusted raster
      contours <- terra::as.polygons(adjusted_altitude_raster, dissolve = TRUE)
      contours <- contours[contours$elevation == 1, ]
      if (nrow(contours) > 0) {
        contours$scientificName <- gsub(" ", "_", species)
        # Export the geographic range as an individual shapefile
      file_name <- paste0(gsub(" ", "_", species), ".shp")
      file_path <- file.path("D:/Doctorado/Tesis/GBIF/contorno", file_name)
      writeVector(contours, file_path, filetype = "ESRI Shapefile", overwrite = TRUE)
      cat("Contour saved as:", file_path, "\n")
    } else {
      cat("No overlapping polygons found for:", species, "\n")
    }
    }
    # Remove temporary objects to reduce memory usage
    rm(each_species, coordinates, points, polygon_sf, polygon_col, dem_contour, adjusted_altitude_raster, contours)
    gc()  # Call garbage collector
  }
  # Apply the geographic range estimation to all species in the dataset
  for (species in species_subset) {
    process_species(species)
  }
  
  # Merge all species range shapefiles into a single dataset
  shapefiles_list <- list.files("D:/Doctorado/Tesis/GBIF/contorno", pattern = "*.shp", full.names = TRUE)
  shapefiles_combined <- lapply(shapefiles_list, st_read) %>% do.call(rbind, .)
  saveRDS(shapefiles_combined, "C:/Users/Dell-PC/Dropbox/CO_DBdata/geographic_range/geographic_range.rds")
  #