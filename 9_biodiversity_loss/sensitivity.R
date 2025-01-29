library(brms)
library(sf)
library(raster)
library(terra)
library(dplyr)
library(doParallel)
library(foreach)
#
# db data set and ecoregions
  db7 <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db7_combined.rds")
  lambda <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/Analysis/get_posterior/posterior_lambda.rds")
  ecoreg <- st_read("D:/Capas/America/ecoregions/ecoreg.shp")
  db_mod_abundance <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod_abundance.rds")
  AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#
# grid (2 km)
  # study_area <- st_transform(ecoreg, crs = AEAstring) 
  # bbox <- st_bbox(study_area)
  # grid_2km <- st_make_grid(study_area, cellsize = c(2000, 2000), square = TRUE)
  # grid_clipped <- st_intersection(grid, study_area)
  # grid_sf <- st_sf(geometry = grid_clipped)
  # grid_sf <- st_join(study_grid, study_area["OBJECTID"], join = st_intersects)
   # st_write(grid_sf, "D:/Capas/America/grid/grid_2km_ID.shp")
  grid_2km <- st_read("D:/Capas/America/grid/grid_2km_ID.shp")
  # grid_centroides <- st_centroid(grid_2km)
  # grid_coords <- st_coordinates(grid_centroides)
   #
   #
  lambda_exp <- lapply(lambda, function(df) {
    df$forest_abundance <- exp(df$log_lambda_0)
    df$pasture_abundance <- exp(df$log_lambda_pasture_offset)
    df <- df[, !(colnames(df) %in% c("log_lambda_0", "log_lambda_pasture_offset"))]
    return(df)
  })
#
# new data for predictions
  unique_subregions <- unique(db7$subregion_species)
  unique_clusters <- unique(db7$cluster_species)
  unique_scientificName <- unique(db7$scientificName)
  subregion_repeat <- sample(unique_subregions, length(grid_coords[,1]), replace = TRUE)
  cluster_repeat <- sample(unique_clusters, length(grid_coords[,1]), replace = TRUE)
  scientificName_repeat <- sample(unique_scientificName, length(grid_coords[,1]), replace = TRUE)
#
# new data for forest scenario
  grid_data_forest <- data.frame(
    longitude = grid_coords[,1],
    latitude = grid_coords[,2],
    scientificName = scientificName_repeat,
    pasture = 0,  # Pasture, variable binaria
    elev_standard = runif(length(grid_coords[,1]), min = min(db7$elev_standard), max = max(db7$elev_standard)),  # Elevación estándar aleatoria
    elev_standard_squared = runif(length(grid_coords[,1]), min = min(db7$elev_standard_squared), max = max(db7$elev_standard_squared)),  # Elevación estándar al cuadrado
    nest_guild = sample(c("Paracoprid", "Telecoprid"), length(grid_coords[,1]), replace = TRUE),  # Niveles de nest_guild
    diet_range = sample(c("coprophagous", "generalist"), length(grid_coords[,1]), replace = TRUE),  # Dieta
    activity = sample(c("diurnal", "nocturnal", "mixed"), length(grid_coords[,1]), replace = TRUE),  # Actividad
    bodysize = runif(length(grid_coords[,1]), min = min(db7$bodysize), max = max(db7$bodysize)),  # Tamaño del cuerpo
    legratio = runif(length(grid_coords[,1]), min = min(db7$legratio), max = max(db7$legratio)),  # Relación de piernas
    subregion_species = subregion_repeat,  # Asignación aleatoria de subregiones
    cluster_species = cluster_repeat  # Asignación aleatoria de clusters
   )
#
# new data for pature scenario  
  grid_data_pasture <- data.frame(
    longitude = grid_coords[,1],
    latitude = grid_coords[,2],
    scientificName = scientificName_repeat,
    pasture = 1,  # Pasture, variable binaria
    elev_standard = runif(length(grid_coords[,1]), min = min(db7$elev_standard), max = max(db7$elev_standard)),  # Elevación estándar aleatoria
    elev_standard_squared = runif(length(grid_coords[,1]), min = min(db7$elev_standard_squared), max = max(db7$elev_standard_squared)),  # Elevación estándar al cuadrado
    nest_guild = sample(c("Paracoprid", "Telecoprid"), length(grid_coords[,1]), replace = TRUE),  # Niveles de nest_guild
    diet_range = sample(c("coprophagous", "generalist"), length(grid_coords[,1]), replace = TRUE),  # Dieta
    activity = sample(c("diurnal", "nocturnal", "mixed"), length(grid_coords[,1]), replace = TRUE),  # Actividad
    bodysize = runif(length(grid_coords[,1]), min = min(db7$bodysize), max = max(db7$bodysize)),  # Tamaño del cuerpo
    legratio = runif(length(grid_coords[,1]), min = min(db7$legratio), max = max(db7$legratio)),  # Relación de piernas
    subregion_species = subregion_repeat,  # Asignación aleatoria de subregiones
    cluster_species = cluster_repeat  # Asignación aleatoria de clusters
  )
  
# Precitions forest scenario for each cell in the grid 
  num_cores <- detectCores() - 2  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  chunk_size <- 3000
  grid_chunks <- split(grid_data_forest, ceiling(seq_along(1:nrow(grid_data_forest)) / chunk_size))
  predictions_forest <- foreach(chunk = grid_chunks, .combine = 'c', .packages = c('brms')) %dopar% {
  pred_chunk <- posterior_predict(db_mod_abundance, newdata = chunk, re_formula=NULL)
  write.table(pred_chunk, file = "./predictions_forest.csv", append = TRUE, col.names = FALSE, row.names = FALSE)
  gc()
  }
  stopCluster(cl)
# prediction pasture scenario for each cell in the grid 
  num_cores <- detectCores() - 2  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  chunk_size <- 3000
  grid_chunks <- split(grid_data_pasture, ceiling(seq_along(1:nrow(grid_data_pasture)) / chunk_size))
  predictions_pasture <- foreach(chunk = grid_chunks, .combine = 'c', .packages = c('brms')) %dopar% {
    pred_chunk <- posterior_predict(db_mod_abundance, newdata = chunk, re_formula=NULL)
    write.table(pred_chunk, file = "./predictions_pasture.csv", append = TRUE, col.names = FALSE, row.names = FALSE)
    gc()
  }
  stopCluster(cl)
   
  
  predictions_forest <- read.csv("./predictions_forest.csv", sep=" ")
   saveRDS(predictions_forest, "./forest_pred.rds")
  p <- readRDS("./forest_pred.rds")
   predictions_pasture <- read.csv("./predictions_pasture.csv", sep=" ")
   saveRDS(predictions_forest, "./pasture_pred.rds")
  
  grid_2km_forest <- grid_2km
  grid_2km_forest$pasture <- 0
  grid_2km_forest$forest_pred <- rowMeans(predictions_forest) 
  
  grid_2km_pasture <- grid_2km
  grid_2km_pasture$pasture <- 1
  grid_2km_pasture$forest_pred <- rowMeans(predictions_pasture) 
  

# Cálculo de sensibilidad por especie
# Para cada pixel, calcula la relación de abundancia entre bosque 
# y potrero. Calcular tasa bosque/potrero
ratio <- predictions %>%
  group_by(species, location) %>%
  summarise(
    abundance_forest = mean(prediction[pasture == "forest"]),
    abundance_pasture = mean(prediction[pasture == "pasture"]),
    sensitivity = abundance_forest / abundance_pasture
  )

# Cálculo de percentiles (25, 50, 75)
# A nivel de cada subregión y escala global, calcula los percentiles de 
# sensibilidad por especie.
percentiles <- ratio %>%
  group_by(species, region) %>%
  summarise(
    low_sensitivity = quantile(sensitivity, 0.25),
    median_sensitivity = quantile(sensitivity, 0.50),
    high_sensitivity = quantile(sensitivity, 0.75)
  )

# Agregación por subregiones biogeográficas
# Agrega las predicciones de biodiversidad a las 13 regiones biogeográficas y a
# todo el área de estudio.
  # Agregar predicciones por regiones biogeográficas
  regional_sensitivity <- ratio %>%
    group_by(region) %>%
    summarise(mean_sensitivity = mean(sensitivity))
  
  # Sensibilidad acumulada al agregar subregiones
  cumulative_sensitivity <- regional_sensitivity %>%
    arrange(region) %>%
    mutate(cumulative_mean = cummean(mean_sensitivity))

# Estimación de biodiversidad relativa a escala pan-Colombia
# Compara la sensibilidad promedio al agregar diferentes cantidades de regiones 
# biogeográficas con el promedio pan-Colombiano.
  pan_colombia_mean <- mean(regional_sensitivity$mean_sensitivity)
  sensitivity_relative <- cumulative_sensitivity %>%
    mutate(relative_to_pan = (cumulative_mean - pan_colombia_mean) / pan_colombia_mean * 100)

# Visualizacion
#  * Percentiles de sensibilidad por especie en diferentes escalas.
#  * Cambios acumulativos de sensibilidad al agregar regiones biogeográficas.
  
  library(ggplot2)
  
  # Percentiles de sensibilidad
  ggplot(percentiles, aes(x = region, y = median_sensitivity)) +
    geom_line() +
    geom_ribbon(aes(ymin = low_sensitivity, ymax = high_sensitivity), alpha = 0.2) +
    labs(title = "Sensitivity percentiles across regions",
         x = "Region", y = "Sensitivity") +
    theme_minimal()
  
  # Sensibilidad acumulativa
  ggplot(sensitivity_relative, aes(x = region, y = relative_to_pan)) +
    geom_line() +
    labs(title = "Relative sensitivity compared to pan-Colombia",
         x = "Number of regions", y = "Relative sensitivity (%)") +
    theme_minimal()