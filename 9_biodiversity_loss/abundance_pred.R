library(brms)
library(sf)
library(raster)
library(terra)
library(dplyr)

db7 <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db7_combined.rds")
lambda <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/Analysis/get_posterior/posterior_lambda.rds")
ecoreg <- st_read("D:/Capas/America/ecoregions/ecoreg.shp")
db_mod_abundance <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod_abundance.rds")
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#
# posterior predict abuncance coefficients based of draws
lambda_exp <- lapply(lambda, function(df) {
  df$forest_abundance <- exp(df$log_lambda_0)
  df$pasture_abundance <- exp(df$log_lambda_pasture_offset)
  df <- df[, !(colnames(df) %in% c("log_lambda_0", "log_lambda_pasture_offset"))]
  return(df)
})

# grid (2 km)
# study_area <- st_transform(ecoreg, crs = AEAstring) 
# bbox <- st_bbox(study_area)
# grid_2km <- st_make_grid(study_area, cellsize = c(2000, 2000), square = TRUE)
# grid_clipped <- st_intersection(grid, study_area)
# grid_sf <- st_sf(geometry = grid_clipped)
# grid_sf <- st_join(study_grid, study_area["OBJECTID"], join = st_intersects)
# st_write(grid_sf, "D:/Capas/America/grid/grid_2km_ID.shp")
# grid_2km <- st_set_crs(grid_2km, AEAstring)
# grid_2km <- st_transform(grid_2km, crs = 4326)
# st_write(grid_2km, "D:/Capas/America/grid/grid_2km_ID_WGS84.shp")

grid_2km <- st_read("D:/Capas/America/grid/grid_2km_ID.shp")
grid_2km$unique_id <- seq_len(nrow(grid_2km))
species <- lambda[[1]]$scientificName

# predictions by grid_2km
n_predictions <- 16

# predictions list by forest and pasture
forest_pred <- list()
pasture_pred <- list()

# Iterar sobre las celdas de la cuadrícula (tienes 229856 celdas)
for (i in 1:nrow(grid_2km)) {
  
  # Crear un dataframe vacío para almacenar las 16 predicciones por especie
  cell_forest_pred <- data.frame(scientificName = lambda[[1]]$scientificName,
                                 matrix(NA, ncol = n_predictions, nrow = 243))
  colnames(cell_forest_pred)[2:17] <- paste0("Prediction_", 1:n_predictions)
  
  cell_pasture_pred <- data.frame(scientificName = lambda[[1]]$scientificName,
                                  matrix(NA, ncol = n_predictions, nrow = 243))
  colnames(cell_pasture_pred)[2:17] <- paste0("Prediction_", 1:n_predictions)
  
  # Iterar sobre los draws (3000 draws)
  for (iter in 1:3000) {
    
    # Obtener los valores de log_lambda_0 y log_lambda_pasture_offset para cada especie en esta iteración
    log_lambda_0 <- lambda[[iter]]$log_lambda_0  # Abundancia para bosque
    log_lambda_pasture_offset <- lambda[[iter]]$log_lambda_pasture_offset  # Abundancia para potrero
    
    # Transformar los log_lambda a abundancia (exponencial)
    forest_abundance <- exp(log_lambda_0)  # Abundancia para bosque
    pasture_abundance <- exp(log_lambda_pasture_offset)  # Abundancia para potrero
    
    # Asignar las 16 predicciones para cada especie en esta celda
    for (k in 1:243) {  # 293 especies
      cell_forest_pred[k, (2 + iter)] <- forest_abundance[k]  # Guardar predicción para bosque
      cell_pasture_pred[k, (2 + iter)] <- pasture_abundance[k]  # Guardar predicción para potrero
    }
  }
  
  # Almacenar los resultados en las listas
  forest_pred[[i]] <- cell_forest_pred
  pasture_pred[[i]] <- cell_pasture_pred
}
saveRDS(forest_pred, "C:/Users/Dell-PC/Dropbox/CO_DBdata/Analysis/predictions/forest_pred.rds")
saveRDS(pasture_pred, "C:/Users/Dell-PC/Dropbox/CO_DBdata/Analysis/predictions/pasture_pred.rds")

################################################################################
