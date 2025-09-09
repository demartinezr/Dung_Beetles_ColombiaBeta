setwd("C:/Users/PC/Dropbox/CO_DBdata")

################################################################################
####################### Grids with varying pixel size ##########################
library(dggridR)
library(sf)
library(raster)
library(dplyr)
library(data.table)
library(purrr)
library(ggplot2)
library(ggridges)
library(gridExtra)
setwd("C:/Users/PC/Dropbox/CO_DBdata")

ecoregions_all <- as.data.frame(readRDS("./ecoregions_all.rds"))
source("F:/repositorio/Dung_Beetles_ColombiaBeta/9_biodiversity_loss/db_compute_loss.R")

# processing by draw
ecoregions_tmp <- ecoregions_all %>%
  select(ecoregions, lon, lat, scientificName, pasture, abun__draw_3000) %>%
  mutate(cell_id = row_number()) %>%
  select(cell_id, everything()) 
row.names(ecoregions_tmp) <- NULL

forest_tmp <- ecoregions_tmp %>%
  filter(pasture == 0) %>%
  select(-pasture)

pasture_tmp <- ecoregions_tmp %>%
  filter(pasture == 1) %>%
  select(-pasture)

# Convertir a formato wide
forest_data <- dcast(as.data.table(forest_tmp), 
                     ecoregions + lon + lat ~ scientificName, 
                     value.var = "abun__draw_3000", fill = 0)

pasture_data <- dcast(as.data.table(pasture_tmp), 
                      ecoregions + lon + lat ~ scientificName, 
                      value.var = "abun__draw_3000", fill = 0)

# get local ratios
cell_ratios <- get_avg_cell_ratios(forest_data, pasture_data, cutoff_type="absolute", cutoff=1)

# Assign each point to a hexagonal cell on a grid of variable resolution
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
points_latlong <- st_as_sf(forest_data[,1:3], coords = c("lon", "lat"), crs = 4326)
points_coords <- as.data.frame(st_coordinates(points_latlong))

# Assign each point to a hexagonal cell on a grid of variable resolution
dggs_4 <- dgconstruct(res=4)
dggs_5 <- dgconstruct(res=5)
dggs_6 <- dgconstruct(res=6)
dggs_7 <- dgconstruct(res=7)
dggs_8 <- dgconstruct(res=8)
dggs_9 <- dgconstruct(res=9)
dggs_10 <- dgconstruct(res=10)

cells_4 <- dgGEO_to_SEQNUM(dggs_4, points_coords$X, points_coords$Y)$seqnum
cells_5 <- dgGEO_to_SEQNUM(dggs_5, points_coords$X, points_coords$Y)$seqnum
cells_6 <- dgGEO_to_SEQNUM(dggs_6, points_coords$X, points_coords$Y)$seqnum
cells_7 <- dgGEO_to_SEQNUM(dggs_7, points_coords$X, points_coords$Y)$seqnum
cells_8 <- dgGEO_to_SEQNUM(dggs_8, points_coords$X, points_coords$Y)$seqnum
cells_9 <- dgGEO_to_SEQNUM(dggs_9, points_coords$X, points_coords$Y)$seqnum
cells_10 <- dgGEO_to_SEQNUM(dggs_10, points_coords$X, points_coords$Y)$seqnum

uc_4 <- unique(cells_4)
uc_5 <- unique(cells_5)
uc_6 <- unique(cells_6)
uc_7 <- unique(cells_7)
uc_8 <- unique(cells_8)
uc_9 <- unique(cells_9)
uc_10 <- unique(cells_10)


# Get average (mean) pointwise medians and the regional medians for each grid cell
cell_pointwise_logratios_4 <- cell_logratios_4 <- cell_richness_4 <- cell_pointwise_richness_4 <- cell_numbers_4 <- cell_beta_4 <- vector()
cell_pointwise_logratios_4_expanded <- cell_logratios_4_expanded <- cell_beta_4_expanded <- rep(NA, nrow(points_latlong))
for(i in 1:length(uc_4)){
  print(i)
  grl <- get_regional_ratios(forest_data, pasture_data, cutoff_type="absolute", cutoff=1, cell_positions = which(cells_4 == uc_4[i]))
  cell_logratios_4[i] <- grl$med_logratio
  cell_richness_4[i] <- grl$n
  cell_logratios_4_expanded[which(cells_4 == uc_4[i])] <- cell_logratios_4[i]
  cell_pointwise_logratios_4[i] <- mean(cell_ratios$med_logratio[which(cells_4 == uc_4[i])], na.rm=T)
  cell_pointwise_richness_4[i] <- mean(cell_ratios$n[cells_4 == uc_4[i] & !is.na(cell_ratios$avg_logratio)])
  cell_beta_4[i] <- cell_richness_4[i] / cell_pointwise_richness_4[i]
  cell_beta_4_expanded[which(cells_4 == uc_4[i])] <- cell_beta_4[i]
  cell_pointwise_logratios_4_expanded[which(cells_4 == uc_4[i])] <- cell_pointwise_logratios_4[i]
  cell_numbers_4[i] <- sum(cells_4 == uc_4[i])
}
raster_regional_4 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_logratios_4_expanded))
raster_pointwise_4 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_pointwise_logratios_4_expanded))
raster_beta_4 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_beta_4_expanded))
raster_difference_4 <- raster_regional_4 - raster_pointwise_4


cell_pointwise_logratios_5 <- cell_logratios_5 <- cell_richness_5 <- cell_pointwise_richness_5 <- cell_numbers_5 <- cell_beta_5 <- vector()
cell_pointwise_logratios_5_expanded <- cell_logratios_5_expanded <- cell_beta_5_expanded <- rep(NA, nrow(points_latlong))
for(i in 1:length(uc_5)){
  print(i)
  grl <- get_regional_ratios(forest_data, pasture_data, cutoff_type="absolute", cutoff=1, cell_positions = which(cells_5 == uc_5[i]))
  cell_logratios_5[i] <- grl$med_logratio
  cell_richness_5[i] <- grl$n
  cell_logratios_5_expanded[which(cells_5 == uc_5[i])] <- cell_logratios_5[i]
  cell_pointwise_logratios_5[i] <- mean(cell_ratios$med_logratio[which(cells_5 == uc_5[i])], na.rm=T)
  cell_pointwise_richness_5[i] <- mean(cell_ratios$n[cells_5 == uc_5[i] & !is.na(cell_ratios$avg_logratio)])
  cell_beta_5[i] <- cell_richness_5[i] / cell_pointwise_richness_5[i]
  cell_beta_5_expanded[which(cells_5 == uc_5[i])] <- cell_beta_5[i]
  cell_pointwise_logratios_5_expanded[which(cells_5 == uc_5[i])] <- cell_pointwise_logratios_5[i]
  cell_numbers_5[i] <- sum(cells_5 == uc_5[i])
}
raster_regional_5 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_logratios_5_expanded))
raster_pointwise_5 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_pointwise_logratios_5_expanded))
raster_beta_5 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_beta_5_expanded))
raster_difference_5 <- raster_regional_5 - raster_pointwise_5


cell_pointwise_logratios_6 <- cell_logratios_6 <- cell_richness_6 <- cell_pointwise_richness_6 <- cell_numbers_6 <- cell_beta_6 <-  vector()
cell_pointwise_logratios_6_expanded <- cell_logratios_6_expanded <- cell_beta_6_expanded <- rep(NA, nrow(points_latlong))
for(i in 1:length(uc_6)){
  print(i)
  grl <- get_regional_ratios(forest_data, pasture_data, cutoff_type="absolute", cutoff=1, cell_positions = which(cells_6 == uc_6[i]))
  cell_logratios_6[i] <- grl$med_logratio
  cell_richness_6[i] <- grl$n
  cell_logratios_6_expanded[which(cells_6 == uc_6[i])] <- cell_logratios_6[i]
  cell_pointwise_logratios_6[i] <- mean(cell_ratios$med_logratio[which(cells_6 == uc_6[i])], na.rm=T)
  cell_pointwise_richness_6[i] <- mean(cell_ratios$n[cells_6 == uc_6[i] & !is.na(cell_ratios$avg_logratio)])
  cell_beta_6[i] <- cell_richness_6[i] / cell_pointwise_richness_6[i]
  cell_beta_6_expanded[which(cells_6 == uc_6[i])] <- cell_beta_6[i]
  cell_pointwise_logratios_6_expanded[which(cells_6 == uc_6[i])] <- cell_pointwise_logratios_6[i]
  cell_numbers_6[i] <- sum(cells_6 == uc_6[i])
}
raster_regional_6 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_logratios_6_expanded))
raster_pointwise_6 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_pointwise_logratios_6_expanded))
raster_beta_6 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_beta_6_expanded))
raster_difference_6 <- raster_regional_6 - raster_pointwise_6


cell_pointwise_logratios_7 <- cell_logratios_7 <- cell_richness_7 <- cell_pointwise_richness_7 <- cell_numbers_7 <- cell_beta_7 <-  vector()
cell_pointwise_logratios_7_expanded <- cell_logratios_7_expanded <- cell_beta_7_expanded <- rep(NA, nrow(points_latlong))
for(i in 1:length(uc_7)){
  print(i)
  grl <- get_regional_ratios(forest_data, pasture_data, cutoff_type="absolute", cutoff=1, cell_positions = which(cells_7 == uc_7[i]))
  cell_logratios_7[i] <- grl$med_logratio
  cell_richness_7[i] <- grl$n
  cell_logratios_7_expanded[which(cells_7 == uc_7[i])] <- cell_logratios_7[i]
  cell_pointwise_logratios_7[i] <- mean(cell_ratios$med_logratio[which(cells_7 == uc_7[i])], na.rm=T)
  cell_pointwise_richness_7[i] <- mean(cell_ratios$n[cells_7 == uc_7[i] & !is.na(cell_ratios$avg_logratio)])
  cell_beta_7[i] <- cell_richness_7[i] / cell_pointwise_richness_7[i]
  cell_beta_7_expanded[which(cells_7 == uc_7[i])] <- cell_beta_7[i]
  cell_pointwise_logratios_7_expanded[which(cells_7 == uc_7[i])] <- cell_pointwise_logratios_7[i]
  cell_numbers_7[i] <- sum(cells_7 == uc_7[i])
}
raster_regional_7 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_logratios_7_expanded))
raster_pointwise_7 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_pointwise_logratios_7_expanded))
raster_beta_7 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_beta_7_expanded))
raster_difference_7 <- raster_regional_7 - raster_pointwise_7

cell_pointwise_logratios_8 <- cell_logratios_8 <- cell_richness_8 <- cell_pointwise_richness_8 <- cell_numbers_8 <- cell_beta_8 <-  vector()
cell_pointwise_logratios_8_expanded <- cell_logratios_8_expanded <- cell_beta_8_expanded <- rep(NA, nrow(points_latlong))
for(i in 1:length(uc_8)){
  print(i)
  grl <- get_regional_ratios(forest_data, pasture_data, cutoff_type="relative", cutoff=.2, cell_positions = which(cells_8 == uc_8[i]))
  cell_logratios_8[i] <- grl$med_logratio
  cell_richness_8[i] <- grl$n
  cell_logratios_8_expanded[which(cells_8 == uc_8[i])] <- cell_logratios_8[i]
  cell_pointwise_logratios_8[i] <- mean(cell_ratios$med_logratio[which(cells_8 == uc_8[i])], na.rm=T)
  cell_pointwise_richness_8[i] <- mean(cell_ratios$n[cells_8 == uc_8[i] & !is.na(cell_ratios$avg_logratio)])
  cell_beta_8[i] <- cell_richness_8[i] / cell_pointwise_richness_8[i]
  cell_beta_8_expanded[which(cells_8 == uc_8[i])] <- cell_beta_8[i]
  cell_pointwise_logratios_8_expanded[which(cells_8 == uc_8[i])] <- cell_pointwise_logratios_8[i]
  cell_numbers_8[i] <- sum(cells_8 == uc_8[i])
}
raster_regional_8 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_logratios_8_expanded))
raster_pointwise_8 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_pointwise_logratios_8_expanded))
raster_beta_8 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_beta_8_expanded))
raster_difference_8 <- raster_regional_8 - raster_pointwise_8

cell_pointwise_logratios_9 <- cell_logratios_9 <- cell_richness_9 <- cell_pointwise_richness_9 <- cell_numbers_9 <- cell_beta_9 <-  vector()
cell_pointwise_logratios_9_expanded <- cell_logratios_9_expanded <- cell_beta_9_expanded <- rep(NA, nrow(points_latlong))
for(i in 1:length(uc_9)){
  print(i)
  grl <- get_regional_ratios(forest_data, pasture_data, cutoff_type="relative", cutoff=.2, cell_positions = which(cells_9 == uc_9[i]))
  cell_logratios_9[i] <- grl$med_logratio
  cell_richness_9[i] <- grl$n
  cell_logratios_9_expanded[which(cells_9 == uc_9[i])] <- cell_logratios_9[i]
  cell_pointwise_logratios_9[i] <- mean(cell_ratios$med_logratio[which(cells_9 == uc_9[i])], na.rm=T)
  cell_pointwise_richness_9[i] <- mean(cell_ratios$n[cells_9 == uc_9[i] & !is.na(cell_ratios$avg_logratio)])
  cell_beta_9[i] <- cell_richness_9[i] / cell_pointwise_richness_9[i]
  cell_beta_9_expanded[which(cells_9 == uc_9[i])] <- cell_beta_9[i]
  cell_pointwise_logratios_9_expanded[which(cells_9 == uc_9[i])] <- cell_pointwise_logratios_9[i]
  cell_numbers_9[i] <- sum(cells_9 == uc_9[i])
}
raster_regional_9 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_logratios_9_expanded))
raster_pointwise_9 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_pointwise_logratios_9_expanded))
raster_beta_9 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_beta_9_expanded))
raster_difference_9 <- raster_regional_9 - raster_pointwise_9

cell_pointwise_logratios_10 <- cell_logratios_10 <- cell_richness_10 <- cell_pointwise_richness_10 <- cell_numbers_10 <- cell_beta_10 <-  vector()
cell_pointwise_logratios_10_expanded <- cell_logratios_10_expanded <- cell_beta_10_expanded <- rep(NA, nrow(points_latlong))
for(i in 1:length(uc_10)){
  print(i)
  grl <- get_regional_ratios(forest_data, pasture_data, cutoff_type="relative", cutoff=.2, cell_positions = which(cells_10 == uc_10[i]))
  cell_logratios_10[i] <- grl$med_logratio
  cell_richness_10[i] <- grl$n
  cell_logratios_10_expanded[which(cells_10 == uc_10[i])] <- cell_logratios_10[i]
  cell_pointwise_logratios_10[i] <- mean(cell_ratios$med_logratio[which(cells_10 == uc_10[i])], na.rm=T)
  cell_pointwise_richness_10[i] <- mean(cell_ratios$n[cells_10 == uc_10[i] & !is.na(cell_ratios$avg_logratio)])
  cell_beta_10[i] <- cell_richness_10[i] / cell_pointwise_richness_10[i]
  cell_beta_10_expanded[which(cells_10 == uc_10[i])] <- cell_beta_10[i]
  cell_pointwise_logratios_10_expanded[which(cells_10 == uc_10[i])] <- cell_pointwise_logratios_10[i]
  cell_numbers_10[i] <- sum(cells_10 == uc_10[i])
}
raster_regional_10 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_logratios_10_expanded))
raster_pointwise_10 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_pointwise_logratios_10_expanded))
raster_beta_10 <- rasterFromXYZ(cbind(forest_data[,2:3], cell_beta_10_expanded))
raster_difference_10 <- raster_regional_10 - raster_pointwise_10

raster_to_df <- function(raster, resolution, type) {
  df <- as.data.frame(raster, xy = TRUE)
  
  # Detectar la columna de logratio
  logratio_col <- grep("cell_.*_expanded|layer", colnames(df), value = TRUE)
  
  df %>%
    rename(logratio = all_of(logratio_col)) %>%
    mutate(res = paste0("Resolution_", resolution),
           type = type)
}

# Lista con todos los rasters y su resolución correspondiente
raster_list <- list(
  "4_regional" = raster_regional_4,
  "4_pointwise" = raster_pointwise_4,
  "4_difference" = raster_difference_4,
  "4_beta" = raster_beta_4,  
  
  "5_regional" = raster_regional_5,
  "5_pointwise" = raster_pointwise_5,
  "5_difference" = raster_difference_5,
  "5_beta" = raster_beta_5,  
  
  "6_regional" = raster_regional_6,
  "6_pointwise" = raster_pointwise_6,
  "6_difference" = raster_difference_6,
  "6_beta" = raster_beta_6,  
  
  "7_regional" = raster_regional_7,
  "7_pointwise" = raster_pointwise_7,
  "7_difference" = raster_difference_7,
  "7_beta" = raster_beta_7
  
  #    "8_regional" = raster_regional_8,
  #    "8_pointwise" = raster_pointwise_8,
  #    "8_difference" = raster_difference_8,
  #    "8_beta" = raster_beta_8,
  
  #    "9_regional" = raster_regional_9,
  #    "9_pointwise" = raster_pointwise_9,
  #    "9_difference" = raster_difference_9,
  #    "9_beta" = raster_beta_9,
  
  #    "10_regional" = raster_regional_10,
  #    "10_pointwise" = raster_pointwise_10,
  #    "10_difference" = raster_difference_10,
  #    "10_beta" = raster_beta_10
)

# Convertir todos los rasters en dataframes
df_list <- lapply(names(raster_list), function(name) {
  parts <- strsplit(name, "_")[[1]]
  resolution <- parts[1]  # Extrae la resolución (4, 8, etc.)
  type <- parts[2]        # Extrae el tipo (regional, pointwise, etc.)
  raster_to_df(raster_list[[name]], resolution, type)
})

# Unir todo en un solo dataframe
df_all <- bind_rows(df_list) %>%
  filter(!is.na(logratio)) %>%
  mutate(res = factor(res, levels = paste0("Resolution_", c(4, 5, 6, 7, 8, 9, 10))))

# Crear subconjuntos
df_subset1 <- df_all %>% filter(type %in% c("regional", "pointwise"))
df_subset2 <- df_all %>% filter(type %in% c("difference", "beta"))

library(tidyr)

# Transformar el formato a "wide" para tener columnas separadas
df_subset2_wide <- df_subset2 %>%
  pivot_wider(names_from = type, values_from = logratio, 
              names_prefix = "") %>%
  rename(excess_loss = difference, beta = beta)

# Graficar la relación entre Beta y Excess Loss
ggplot(df_subset2_wide, aes(x = beta, y = excess_loss)) +
  geom_hex(bins = 30) +  # Mapa de densidad hexagonal
  scale_fill_viridis_c(trans = "log") +  # Escala perceptual
  geom_smooth(method = "lm", color = "black", fill = "gray40", se = TRUE) +  # Línea de regresión con IC
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~res, nrow = 2) +  # Facetas por resolución
  theme_bw(base_size = 14) +
  labs(x = "Beta", y = "Excess regional loss", fill = "Density") +
  #    scale_x_continuous(limits=c(0, 20)) +
  theme(legend.position = "none")

