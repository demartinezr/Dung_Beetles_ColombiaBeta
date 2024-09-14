# Function for obtaining geographic ranges based on SIB_col and collections, 
# ecoregions, and altitudinal ranges (alos palsar DEM 30m) for morphospecies.
#
  setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata")
#
# R packages
    library(raster)
    library(dplyr)
    library(ggplot2)
    library(rgdal)
    library(sp)
    library(rgeos)
    library(terra)
    library(foreach)
    library(doParallel)
#
# Base maps
  ecoreg <- readOGR(dsn = "D:/Capas/World/wwf/Ecoregions2017", layer = "eco_col") # País
  ecoreg <- spTransform(ecoreg, CRS("+proj=longlat +datum=WGS84")) # WGS84
  ecoreg <- as(ecoreg, "Spatial")
  dem_ALOS <- "D:/Capas/coldem30/alos_dem_col.tif"
  dem <- raster(dem_ALOS)
#
# GBIF data
  registros <- readRDS("morphoIAvH_SIB_ele.rds")
  proj_crs <- CRS("+proj=utm +zone=18 +datum=WGS84 +units=m +no_defs")  
  registros <- mutate(registros, identificationRemarks = ifelse(identificationRemarks=="Canthidium_sp_13H", "Canthidium sp. 13H", identificationRemarks))
  registros <- mutate(registros, identificationRemarks = ifelse(identificationRemarks=="Canthon_sp_11H", "Canthon sp. 11H", identificationRemarks))
  registros <- mutate(registros, identificationRemarks = ifelse(identificationRemarks=="Canthon_sp_18H", "Canthon sp. 18H", identificationRemarks))
  registros <- mutate(registros, identificationRemarks = ifelse(identificationRemarks=="Canthon_sp_19H", "Canthon sp. 19H", identificationRemarks))
  registros <- mutate(registros, identificationRemarks = ifelse(identificationRemarks=="Canthon_sp_22H", "Canthon sp. 22H", identificationRemarks))
  registros <- mutate(registros, identificationRemarks = ifelse(identificationRemarks=="Deltochilum_sp_18H", "Deltochilum sp. 18H", identificationRemarks))
  registros <- mutate(registros, identificationRemarks = ifelse(identificationRemarks=="Onthophagus_sp_08H", "Onthophagus sp. 08H", identificationRemarks))
  registros <- mutate(registros, identificationRemarks = ifelse(identificationRemarks=="Uroxys_sp_03H", "Uroxys sp. 03H", identificationRemarks))
  registros <- mutate(registros, identificationRemarks = ifelse(identificationRemarks=="Uroxys_sp_06H", "Uroxys sp. 06H", identificationRemarks))
  registros <- mutate(registros, identificationRemarks = ifelse(identificationRemarks=="Uroxys_sp_08H", "Uroxys sp. 08H", identificationRemarks))
#
# Species subset
  G1 <- unique(registros$identificationRemarks)[c(34, 37, 39, 43, 45, 46)] 
#
# Function to process each species
  process_species <- function(especie) 
{
    cada_especie <- registros[registros$identificationRemarks == especie, ]
  #
  # GBIF records by species
    coordenadas <- cada_especie[, c("decimalLongitude", "decimalLatitude")]
    coordenadas <- na.omit(coordenadas)
    coordenadas <- coordenadas[!duplicated(coordenadas), ]
    puntos <- SpatialPoints(coordenadas, proj4string = CRS("+proj=longlat +datum=WGS84"))
    puntos_col <- over(puntos, as(ecoreg, "SpatialPolygons"))
    puntos_col <- !is.na(puntos_col)
    puntos_col <- puntos[puntos_col, ]
    if (length(puntos_col) == 0) {
      cat("No records in Colombia for:", especie, "\n")
      return(NULL)
    }
    coordenadas <- as.data.frame(puntos_col)
  #
  # Polygon based on the number of records 1, 2, 3 ...
    if (nrow(coordenadas) <= 2) {
      if (nrow(coordenadas) == 1) {
        punto <- SpatialPoints(coordenadas, proj4string = CRS("+proj=longlat +datum=WGS84"))
        punto <- spTransform(punto, CRS("+proj=longlat +datum=WGS84 +no_defs"))
        poligono_sp <- gBuffer(punto, width = 0.5)
      } else {
        linea <- SpatialLines(list(Lines(list(Line(coordenadas)), "1")), proj4string = CRS("+proj=longlat +datum=WGS84"))
        linea <- spTransform(linea, CRS("+proj=longlat +datum=WGS84 +no_defs"))
        poligono_sp <- gBuffer(linea, width = 0.5)
        proj4string(poligono_sp) <- CRS("+proj=longlat +datum=WGS84")
      }
    } else {
      indices_poligono <- chull(coordenadas$decimalLongitude, coordenadas$decimalLatitude)
      poligono <- coordenadas[indices_poligono, ]
      poligono_sp <- SpatialPolygons(list(Polygons(list(Polygon(poligono)), "1")))
      proj4string(poligono_sp) <- CRS("+proj=longlat +datum=WGS84")
    }
    pol_sol <- tryCatch(gIntersects(ecoreg, poligono_sp, byid = TRUE), error = function(e) NULL)
    if (!is.null(pol_sol)) {
      indice_sol <- which(pol_sol)
      if (length(indice_sol) > 0) {
        sol_com <- ecoreg[indice_sol, ]
        sol_com <- spTransform(sol_com, proj_crs)
        sol_com <- gBuffer(sol_com, byid = TRUE, width = 0)
        sol_com <- spTransform(sol_com, CRS("+proj=longlat +datum=WGS84 +no_defs"))
        contorno <- gUnaryUnion(sol_com)
  #
  # Adjust contour to the altitudinal limits of the species
    dem_contorno <- mask(dem, contorno)
    altura_max <- max(cada_especie$alos_ele, na.rm = TRUE) + 100
    altura_min <- max(min(cada_especie$alos_ele, na.rm = TRUE) - 100, 0)
    raster_altitud_ajustada <- dem_contorno >= altura_min & dem_contorno <= altura_max
    spat_raster <- terra::rast(raster_altitud_ajustada)
    #     
    # Raster to polygon
      dem_ajustado <- terra::as.polygons(spat_raster)
      poligonos_valor_1 <- dem_ajustado[dem_ajustado$layer == 1, ]
      poligono_union <- terra::aggregate(poligonos_valor_1)
  #
  # Save
    nombre_archivo <- paste0(gsub(" ", "._", especie), ".shp")
    ruta_archivo <- file.path("D:/Doctorado/Tesis/GBIF/contorno", nombre_archivo)
    writeVector(poligono_union, ruta_archivo, filetype = "ESRI Shapefile")
    return(nombre_archivo)
    } else {
      cat("No overlapping polygons found for:", especie, "\n")
      return(NULL)
    }
  } else {
    cat("Error intersecting polygons for:", especie, "\n")
    return(NULL)
  }
}
#
# Set up parallel backend to use 4 cores
  cl <- makeCluster(4)
  registerDoParallel(cl)
  #
  # Parallel processing
   resultados_especies <- foreach(especie = G1, .combine = 'c', .packages = c('raster', 'rgdal', 'sp', 'rgeos', 'terra', 'dplyr')) %dopar% 
   {
   process_species(especie)
   }
  #
  # Stop the cluster
    stopCluster(cl)
  #
  # Convert the list of results to a named list
  names(species_result) <- G1
  #
  # Filter out NULL results
  species_result <- species_result[!sapply(species_result, is.null)]

