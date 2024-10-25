# ALOS-Palsar dem 30m assembly
#
  setwd("D:/Users/DELL/Documents/Capas/america_S_dem")
  library(terra)
  # List all the TIFF files you want to merge
  tif_files <- list.files(pattern = "\\.tif$")
  # Read the TIFF files into a SpatRaster object
  rasters <- lapply(tif_files, rast)
  # Merge all the rasters into a single one
  merged_raster <- do.call(merge, rasters)
  # Save the merged raster to a new TIFF file
  writeRaster(merged_raster, "srtm90m.tif", overwrite = TRUE)
  # Plot the merged raster to visualize it
  plot(merged_raster)
#
# strm dem 90m
    setwd("D:/Capas/America/dem/srtm")
    # Load the required library
    library(terra)
    # List all the TIFF files you want to merge
    tif_files <- list.files(pattern = "\\.tif$")
    # Read the TIFF files into a SpatRaster object
    rasters <- lapply(tif_files, rast)
    # Merge all the rasters into a single one
    merged_raster <- do.call(merge, rasters)
    # Save the merged raster to a new TIFF file
    writeRaster(merged_raster, "srtm90m.tif", overwrite = TRUE)
        # Plot the merged raster to visualize it
    plot(merged_raster)


  
  