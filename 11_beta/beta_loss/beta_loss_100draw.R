################################################################################
############# hexagons approach: Grids with varying pixel size #################
setwd("C:/Users/PC/Dropbox/CO_DBdata")
# Spatial and raster processing
library(dggridR)
library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(patchwork)

# Compute local and regional diversity loss across multiple spatial resolutions and posterior draws

# Assign each point to a hexagonal cell on a grid of variable resolution
# AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
# forest_data <- ecoregions_all %>% filter(pasture==1)
# points_latlong <- forest_data %>%  select(ecoregions, lat, lon) %>%  distinct() %>%  st_as_sf(coords = c("lon", "lat"), crs = 4326)
# points_coords <- as.data.frame(st_coordinates(points_latlong))
# saveRDS(points_coords, "./diversity_loss/points_coords.rds")

# Load custom functions from the 'compute_loss.R' script
source("F:/repositorio/Dung_Beetles_ColombiaBeta/9_biodiversity_loss/db_compute_loss.R")

# define resolutions, draws and data
resolutions <- c(5:10)
#  draws <- seq(300, 3000, by = 300)
draws <-  30 * c(1:100) 
ecoregions_all <- as.data.frame(readRDS("./ecoregions_100.rds"))
points_coords <- readRDS("./diversity_loss/points_coords.rds")

# Generate the forest and pasture lists in wide format for each draw
# This function returns a list containing two elements: forest_list and pasture_list
draw_data <- generate_forest_pasture_lists(ecoregions_all, draws)
#  saveRDS(draw_data, "./diversity_loss/draw_data_100.rds")
draw_data <- readRDS("./diversity_loss/draw_data_100.rds")
# Assign the results to global variables for further analysis
forest_list <- draw_data$forest_list
pasture_list <- draw_data$pasture_list

# Build DGGS grids for each resolution
dggs <- lapply(resolutions, function(x) dgconstruct(res = x))
names(dggs) <- paste0("dggs_", resolutions)

# Convert coordinates to cell sequence numbers
cells <- lapply(dggs, function(dg) dgGEO_to_SEQNUM(dg, points_coords$X, points_coords$Y)$seqnum)
uc <- lapply(cells, unique)

# Initialize the list to store raster results for each draw and resolution   
resultados_rasters <- list()

# Main loop through draws and resolutions
for (draw in draws) {
  draw_name <- paste0("draw_", draw)
  resultados_rasters[[draw_name]] <- list()
  
  for (r in seq_along(resolutions)) {
    cat("draw:", draw, "res:", resolutions[r], "\n")
    
    # Run the analysis for each combination of draw and resolution
    raster_layers <- process_draw_resolution(
      draw = draw,
      res_index = r,
      resolutions = resolutions,
      dggs_list = dggs,
      cell_list = cells,
      unique_cells_list = uc,
      forest_list = forest_list,
      pasture_list = pasture_list,
      coords = points_coords)
    
    # Save the results for this resolution
    res_name <- paste0("res_", resolutions[r])
    resultados_rasters[[draw_name]][[res_name]] <- raster_layers
  }
}

saveRDS(resultados_rasters, "./diversity_loss/resultados_rasters_100.rds")

resultados_rasters <- readRDS("./diversity_loss/resultados_rasters_10.rds")
# Extraer datos de todos los draws y resoluciones
df_list <- list()

for (draw in names(resultados_rasters)) {
  for (res in names(resultados_rasters[[draw]])) {
    #        df_list[[paste0(draw, "_", res, "_regional")]] <-
    #          raster_to_df(resultados_rasters[[draw]][[res]]$raster_regional, resolution = res, draw = draw, type = "Regional", colname = "logratios")
    
    #        df_list[[paste0(draw, "_", res, "_pointwise")]] <-
    #          raster_to_df(resultados_rasters[[draw]][[res]]$raster_pointwise, resolution = res, draw = draw, type = "Pointwise", colname = "pointwise")
    
    df_list[[paste0(draw, "_", res, "_difference")]] <-
      raster_to_df(resultados_rasters[[draw]][[res]]$raster_difference, resolution = res, draw = draw, type = "Difference", colname = "layer")
    
    df_list[[paste0(draw, "_", res, "_beta")]] <-
      raster_to_df(resultados_rasters[[draw]][[res]]$raster_beta, resolution = res, draw = draw, type = "Beta", colname = "beta")
  }
}
#  saveRDS(df_list, "./diversity_loss/df_list.rds")  

# from list to dataframe
df_list <- readRDS("./diversity_loss/df_list_DB.rds")

df_all <- bind_rows(df_list) %>%
  filter(!is.na(logratio))
saveRDS(df_list, "./diversity_loss/df_all_DB.rds")  

#    df_subset1 <- df_all %>% filter(type %in% c("Regional", "Pointwise"))
#    df_subset2 <- df_all %>% filter(type %in% c("Difference", "Beta"))

#    ggplot(df_subset1, aes(x = x, y = y, fill = logratio)) +
#      geom_raster() +
#      scale_fill_viridis_c() +
#      facet_grid(type ~ res) + 
#      theme_bw() +
#      labs(title = "Comparison of regional and pointwise Log ratios", 
#          fill = "Log ratio")


ggplot(df_all, aes(x = x, y = y, fill = logratio)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(type ~ res, scales = "free") + 
  theme_bw() +
  labs(title = NULL, 
       fill = "ratio")
################################
# fig 3 

resolutions <- c(5, 7, 10)
dggs <- lapply(resolutions, function(x) dgconstruct(res = x))
names(dggs) <- paste0("dggs_", resolutions)

draw_300 <- resultados_rasters$draw_300

raster_res_5_diff <- draw_300$res_5$raster_difference
raster_res_7_diff <- draw_300$res_7$raster_difference
raster_res_10_diff <- draw_300$res_10$raster_difference

raster_res_5_beta <- draw_300$res_5$raster_beta
raster_res_7_beta <- draw_300$res_7$raster_beta
raster_res_10_beta <- draw_300$res_10$raster_beta


# Función para convertir raster a valores hexagonales para cada resolución
convert_raster_to_sfhex <- function(raster_obj, dggs_obj) {
  # Extraer coordenadas de las celdas del raster
  coords <- xyFromCell(raster_obj, 1:ncell(raster_obj))
  values <- values(raster_obj)
  
  # Filtrar celdas con NA
  valid_idx <- which(!is.na(values))
  coords <- coords[valid_idx, ]
  values <- values[valid_idx]
  
  # Asignar cada punto a un hexágono del sistema DGGS
  hex_ids <- dgGEO_to_SEQNUM(dggs_obj, coords[, 1], coords[, 2])$seqnum
  
  # Promediar valores por hexágono
  df <- data.frame(hex_id = hex_ids, value = values) %>%
    group_by(hex_id) %>%
    summarize(value = mean(value, na.rm = TRUE))
  
  # Obtener la geometría de cada hexágono
  hex_centroids <- dgSEQNUM_to_GEO(dggs_obj, df$hex_id)
  hex_polys <- dgcellstogrid(dggs_obj, df$hex_id)
  
  hex_sf <- st_as_sf(hex_polys)
  hex_sf$value <- df$value
  
  return(hex_sf)
}

# Para Excess Loss
hex_res5_diff <- convert_raster_to_sfhex(raster_res_5_diff, dggs$dggs_5)
hex_res7_diff <- convert_raster_to_sfhex(raster_res_7_diff, dggs$dggs_7)
hex_res10_diff <- convert_raster_to_sfhex(raster_res_10_diff, dggs$dggs_10)

# Para Beta-diversity
hex_res5_beta <- convert_raster_to_sfhex(raster_res_5_beta, dggs$dggs_5)
hex_res7_beta <- convert_raster_to_sfhex(raster_res_7_beta, dggs$dggs_7)
hex_res10_beta <- convert_raster_to_sfhex(raster_res_10_beta, dggs$dggs_10)

colombia <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(name == "Colombia")
colombia <- st_transform(colombia, crs = st_crs(hex_res5_diff))


# Recortar los hexágonos para que solo queden dentro del contorno de Colombia
hex_res5_diff_clipped <- st_intersection(hex_res5_diff, colombia)
hex_res7_diff_clipped <- st_intersection(hex_res7_diff, colombia)
hex_res10_diff_clipped <- st_intersection(hex_res10_diff, colombia)

hex_res5_beta_clipped <- st_intersection(hex_res5_beta, colombia)
hex_res7_beta_clipped <- st_intersection(hex_res7_beta, colombia)
hex_res10_beta_clipped <- st_intersection(hex_res10_beta, colombia)

###### plot ###

p_diff_res5 <- ggplot() +
  geom_sf(data = colombia, fill = "gray", color = "black") +  # Fondo de Colombia
  geom_sf(data = hex_res5_diff_clipped, aes(fill = value), color = NA) +
  ggtitle("(a)") +
  theme_bw() +
  scale_fill_viridis_c(name = "Excess Loss", option = "C", na.value = NA) +
  scale_x_continuous(breaks = seq(-80, -65, length.out = 4))  # Dividir eje X en 4 puntos

p_diff_res7 <- ggplot() +
  geom_sf(data = colombia, fill = "gray", color = "black") +  # Fondo de Colombia
  geom_sf(data = hex_res7_diff_clipped, aes(fill = value), color = NA) +
  ggtitle("(b)") +
  theme_bw() +
  scale_fill_viridis_c(name = "Excess Loss", option = "C", na.value = NA) +
  scale_x_continuous(breaks = seq(-80, -65, length.out = 4))  # Dividir eje X en 4 puntos

p_diff_res10 <- ggplot() +
  geom_sf(data = colombia, fill = "gray", color = "black") +  # Fondo de Colombia
  geom_sf(data = hex_res10_diff_clipped, aes(fill = value), color = NA) +
  ggtitle("(c)") +
  theme_bw() +
  scale_fill_viridis_c(name = "Excess Loss", option = "C", na.value = NA) +
  scale_x_continuous(breaks = seq(-80, -65, length.out = 4))  # Dividir eje X en 4 puntos

# Paneles de Beta-diversity (d–f)
p_beta_res5 <- ggplot() +
  geom_sf(data = colombia, fill = "gray", color = "black") +  # Fondo de Colombia
  geom_sf(data = hex_res5_beta_clipped, aes(fill = value), color = NA) +
  ggtitle("(d)") +
  theme_bw() +
  scale_fill_viridis_c(name = "Beta", option = "D", na.value = NA) +
  scale_x_continuous(breaks = seq(-80, -65, length.out = 4))  # Dividir eje X en 4 puntos

p_beta_res7 <- ggplot() +
  geom_sf(data = colombia, fill = "gray", color = "black") +  # Fondo de Colombia
  geom_sf(data = hex_res7_beta_clipped, aes(fill = value), color = NA) +
  ggtitle("(e)") +
  theme_bw() +
  scale_fill_viridis_c(name = "Beta", option = "D", na.value = NA) +
  scale_x_continuous(breaks = seq(-80, -65, length.out = 4))  # Dividir eje X en 4 puntos

p_beta_res10 <- ggplot() +
  geom_sf(data = colombia, fill = "gray", color = "black") +  # Fondo de Colombia
  geom_sf(data = hex_res10_beta_clipped, aes(fill = value), color = NA) +
  ggtitle("(f)") +
  theme_bw() +
  scale_fill_viridis_c(name = "Beta", option = "D", na.value = NA) +
  scale_x_continuous(breaks = seq(-80, -65, length.out = 4))  # Dividir eje X en 4 puntos

# Crear las filas y la figura final
fila1 <- p_diff_res5 + p_diff_res7 + p_diff_res10
fila2 <- p_beta_res5 + p_beta_res7 + p_beta_res10
fig_4 <- fila1 / fila2

figura_final

ggsave("./fig_3.jpeg", plot = fig_3, width = 8.5, height = 5.5, units = "in",        
       dpi = 300, device = "jpeg")
####################################################
# fig 4
# Transformar el formato a "wide" para tener columnas separadas
df_subset2_wide <- df_all %>%
  pivot_wider(names_from = type, values_from = logratio, 
              names_prefix = "") %>%
  rename(excess_loss = Difference, Beta = Beta)

#    saveRDS(df_subset2_wide, "./diversity_loss/df_subset2_wide.rds")

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

df_subset2_wide <- readRDS("./diversity_loss/df_subset2_wide.rds")

df_subset2_wide <- df_subset2_wide %>%
  mutate(res = case_when(
    res == "res_5"  ~ "70000~km^2",
    res == "res_6"  ~ "23000~km^2",
    res == "res_7"  ~ "7800~km^2",
    res == "res_8"  ~ "2600~km^2",
    res == "res_9"  ~ "860~km^2",
    res == "res_10" ~ "290~km^2",
    TRUE            ~ res
  ))

df_subset2_wide$res <- factor(df_subset2_wide$res, 
                              levels = c("70000~km^2", "23000~km^2", "7800~km^2", 
                                         "2600~km^2", "860~km^2", "290~km^2"))

# 1. Extraer valores únicos de beta observados (no predicción)
beta_vals <- sort(unique(df_subset2_wide$Beta))

# 2. Obtener coeficientes por draw y res
coef_df <- df_subset2_wide %>%
  group_by(res, draw) %>%
  group_split() %>%
  map_dfr(function(df) {
    fit <- lm(excess_loss ~ Beta, data = df)
    tibble(
      res = unique(df$res),
      draw = unique(df$draw),
      intercept = coef(fit)[1],
      slope = coef(fit)[2]
    )
  })

#   saveRDS(coef_df, "./diversity_loss/coef_df.rds")

# 3. Construir rectas por cada beta observado
line_df <- expand_grid(beta = beta_vals, coef_df) %>%
  mutate(pred = intercept + slope * beta)
#    saveRDS(line_df, "./diversity_loss/line_df.rds")

# 4. Promediar las rectas por res y beta
summary_df <- line_df %>%
  group_by(res, beta) %>%
  summarise(
    ymed = mean(pred),
    ymin = quantile(pred, 0.025),
    ymax = quantile(pred, 0.975),
    .groups = "drop"
  )
#    saveRDS(summary_df, "./diversity_loss/summary_df.rds")     
summary_df <- readRDS("./diversity_loss/summary_df.rds")

summary_df <- summary_df %>%
  mutate(res = case_when(
    res == "res_5"  ~ "70000~km^2",
    res == "res_6"  ~ "23000~km^2",
    res == "res_7"  ~ "7800~km^2",
    res == "res_8"  ~ "2600~km^2",
    res == "res_9"  ~ "860~km^2",
    res == "res_10" ~ "290~km^2",
    TRUE            ~ res
  ))

summary_df$res <- factor(summary_df$res, 
                         levels = c("70000~km^2", "23000~km^2", "7800~km^2", 
                                    "2600~km^2", "860~km^2", "290~km^2"))    

# 5. Graficar
fig_4 <- ggplot(df_subset2_wide, aes(x = Beta, y = excess_loss)) +
  geom_hex(bins = 30) +
  scale_fill_viridis_c(trans = "log") +
  geom_ribbon(data = summary_df, aes(x = beta, ymin = ymin, ymax = ymax),
              fill = "gray40", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = summary_df, aes(x = beta, y = ymed),
            color = "black", linewidth = 1, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  facet_wrap(~res, ncol = 3, labeller = label_parsed) +
  theme_bw(base_size = 14) +
  labs(x = "Beta", y = "Excess regional loss", fill = "Density") +
  theme(legend.position = "none")

ggsave("./fig_4.jpeg", plot = fig_4, width = 8.5, height = 5.5, units = "in",        
       dpi = 300, device = "jpeg")
