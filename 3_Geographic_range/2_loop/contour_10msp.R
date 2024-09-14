# Loop for obtaining geographic ranges from GBIF records, ecoregions, and 
# altitudinal ranges (DEM)
setwd("D:/Doctorado/Tesis/GBIF")
#
# R packages
library(raster)
library(dplyr)
library(ggplot2)
library(rgdal)
library(sp)
library(rgeos)
library(terra)
#
# Base maps
  ecoreg <- readOGR(dsn = "D:/Capas/World/wwf/Ecoregions2017", layer = "eco_col") #Pais
  ecoreg <- spTransform(ecoreg, CRS("+proj=longlat +datum=WGS84")) # WGS84
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
# loop
  resultados_especies <- list() # result list
  G1 <- unique(registros$identificationRemarks)[17:20] # species subset
  for (especie in G1) {
  cada_especie <- registros[registros$identificationRemarks == especie, ]
  # GBIF records by specie
    coordenadas <- cada_especie[, c("decimalLongitude", "decimalLatitude")]
    coordenadas <- na.omit(coordenadas)
    coordenadas <- coordenadas[!duplicated(coordenadas), ]
    puntos <- SpatialPoints(coordenadas, proj4string = CRS("+proj=longlat +datum=WGS84"))  # coordinates to SpatialPoints
    puntos_col <- over(puntos, as(ecoreg, "SpatialPolygons")) # Filter points in Colombia
    puntos_col <- !is.na(puntos_col)
    puntos_col <- puntos[puntos_col, ]
    if (length(puntos_col) == 0) {
    cat("No records in Colombia for:", especie, "\n")
    next
    }
    coordenadas <- as.data.frame(puntos_col)
  # polygon based on the number of records 1, 2, 3 ...
    if (nrow(coordenadas) <= 2) {
    if (nrow(coordenadas) == 1) {
      # Buffer on a single coordinate
        punto <- SpatialPoints(coordenadas, proj4string = CRS("+proj=longlat +datum=WGS84"))
        punto <- spTransform(punto, CRS("+proj=longlat +datum=WGS84 +no_defs"))
        poligono_sp <- gBuffer(punto, width = 0.5) # Buffer de 0.01 grades (~1 km)
        } else {
      # Polygon based on two coordinates
        linea <- SpatialLines(list(Lines(list(Line(coordenadas)), "1")), proj4string = CRS("+proj=longlat +datum=WGS84"))
        linea <- spTransform(linea, CRS("+proj=longlat +datum=WGS84 +no_defs"))
        poligono_sp <- gBuffer(linea, width = 0.5) # Buffer de 0.01 grades (~1 km)
        proj4string(poligono_sp) <- CRS("+proj=longlat +datum=WGS84 +units=m +no_defs")
        }
        } else {
      # Polygon based on 3 or more cordinates
        indices_poligono <- chull(coordenadas$decimalLongitude, coordenadas$decimalLatitude)
        poligono <- coordenadas[indices_poligono, ]
        poligono_sp <- SpatialPolygons(list(Polygons(list(Polygon(poligono)), "1")))
        proj4string(poligono_sp) <- CRS("+proj=longlat +datum=WGS84")
        }
        pol_sol <- tryCatch(gIntersects(ecoreg, poligono_sp, byid = TRUE), error = function(e) NULL) # Polygon intersection
        if (!is.null(pol_sol)) {
        indice_sol <- which(pol_sol) # Indices of overlapping polygons
        if (length(indice_sol) > 0) {
        sol_com <- ecoreg[indice_sol, ] # Extract the complete overlapping polygons
        sol_com <- spTransform(sol_com, proj_crs)
        sol_com <- gBuffer(sol_com, byid=TRUE, width=0)
        sol_com <- spTransform(sol_com, CRS("+proj=longlat +datum=WGS84 +no_defs"))
        contorno <- gUnaryUnion(sol_com)
  # Adjust contour to the altitudinal limits of the species
      # Clipping DEM by contours
        dem_contorno <- mask(dem, contorno)
      # Altitudinal limits by species
        altura_max <- max(cada_especie$alos_ele, na.rm = TRUE) + 100
        altura_min <- max(min(cada_especie$alos_ele, na.rm = TRUE) - 100, 0)
        raster_altitud_ajustada <- dem_contorno >= altura_min & dem_contorno <= altura_max # Logical raster marking TRUE for cells within the altitudinal range
        spat_raster <- terra::rast(raster_altitud_ajustada)
      # Raster to polygon
        dem_ajustado <- terra::as.polygons(spat_raster) # SpatRaster to polygons
        poligonos_valor_1 <- dem_ajustado[dem_ajustado$layer == 1, ] # Select only polygons with a value of 1
        poligono_union <- terra::aggregate(poligonos_valor_1) # Merge polygons with a value of 1 into a single polygon
  # Save
    nombre_archivo <- paste0("poligono_", gsub(" ", "_", especie), ".shp")
    ruta_archivo <- file.path("D:/Doctorado/Tesis/GBIF/contorno", nombre_archivo)
    writeVector(poligono_union, ruta_archivo, filetype = "ESRI Shapefile")
  # Save the filename in the species results list
    resultados_especies[[especie]] <- nombre_archivo
    } else {
    cat("No overlapping polygons found for:", especie, "\n")
    }
    } else {
    cat("Error intersecting polygons for:", especie, "\n")
   }
  }

plot(poligono_union)
################# Visualización ####################
  ### Ruta al directorio que contiene los archivos shapefile
  directorio_shapes <- "D:/Doctorado/Tesis/GBIF/contorno"
  ### Obtener la lista de archivos de shapefile en el directorio
  archivos_shapes <- list.files(directorio_shapes, pattern = "\\.shp$", full.names = TRUE)
  ### Crear una lista para almacenar los objetos de los shapefiles
  lista_shapefiles <- list()
  ### Iterar sobre los archivos de shapefile y cargarlos en la lista
  for (archivo_shape in archivos_shapes)
    {
    ### Obtener el nombre de la especie del archivo shape
    nombre_especie <- tools::file_path_sans_ext(basename(archivo_shape))
    ### Cargar el shapefile y agregarlo a la lista
    shapefile <- readOGR(dsn = archivo_shape, layer = tools::file_path_sans_ext(basename(archivo_shape)))
    lista_shapefiles[[nombre_especie]] <- shapefile
    }

  sp <- data.frame(subset(registros, identificationRemarks == "Deltochilum_sp_18H"))
  coordenadas <- sp[,c("decimalLongitude", "decimalLatitude")]
  coordenadas <- na.omit(coordenadas)
  coordenadas <- coordenadas[!duplicated(coordenadas), ]
par(mar=c(2,2,1,1)) 
plot(lista_shapefiles[[15]], main = names(lista_shapefiles)[1])
plot(ecoreg)
plot(puntos_en_colombia, add = TRUE, pch = 2)
points(coordenadas$decimalLongitude, coordenadas$decimalLatitude, col = "blue", pch = 20, cex=2)


