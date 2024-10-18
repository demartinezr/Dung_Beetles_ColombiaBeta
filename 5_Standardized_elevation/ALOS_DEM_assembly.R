setwd("D:/Users/DELL/Documents/Capas/america_S_dem")

library(raster)

directorio <- "D:/Users/DELL/Documents/Capas/america_S_dem/"

# Encuentra todos los archivos TIFF en el directorio especificado
archivos_tif <- list.files(path = directorio, pattern = "\\.tif$", full.names = TRUE)

# Comprueba si se han encontrado los archivos TIFF
print(archivos_tif)

# Verifica que haya al menos un archivo encontrado
if (length(archivos_tif) > 0) {
  # Carga todos los archivos TIFF en una lista
  raster_list <- lapply(archivos_tif, raster)
  
  # Fusiona los objetos raster en uno solo
  dem_combined <- merge(raster_list)
  
  # Guarda el DEM combinado en un archivo TIFF
  writeRaster(dem_combined, filename = "dem_combined.tif", format = "GTiff")
} else {
  print("No se encontraron archivos TIFF en el directorio especificado.")
  }
  
  
  