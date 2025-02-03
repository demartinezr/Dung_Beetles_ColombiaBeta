setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata")

library(sf)
library(dplyr)
library(ggplot2)

# Abundance predictions for 243 species of dung beetles with 10 iterations in 
# pasture and forest
db_predictions <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/species_predictions.rds")

# Study area ecoregions and grid
study_area <- st_read("D:/Capas/America/ecoregions/ecoreg.shp")
study_area <- st_make_valid(study_area)

# Create a list to store sf dataframes by ecoregion
ec_points <- vector("list", length = nrow(study_area))
names(ec_points) <- study_area$ECO_NAME  # Use ecoregion names

# Iterate over each ecoregion
for (i in seq_len(nrow(study_area))) {
  polygon <- study_area[i, ]
  # Filter points of all species that fall within the ecoregion
  points_in_ecoregion <- do.call(rbind, lapply(db_predictions, function(df) {
    st_filter(df, polygon)
  }))
  # Store the result in the list
  ec_points[[i]] <- points_in_ecoregion
}
saveRDS(ec_points, "./Analysis/mean_abundance/ecoregions_predictions.rds")
gc()
################################
#
# comparisons
eco_predictions <- readRDS("./Analysis/mean_abundance/ecoregions_predictions.rds")
sta_marta <- eco_predictions[["Santa Marta montane forests"]]
magdalena_dry <- eco_predictions[["Magdalena Valley dry forests"]]
cauca_montane <- eco_predictions[["Cauca Valley montane forests"]]
EC_montane <- eco_predictions[["Eastern Cordillera Real montane forests"]]
villavicencio_dry <- eco_predictions[["Apure-Villavicencio dry forests"]]

# function to calculate the multiplicative change of abundance by species in 
# 10 posterior iterations

calcular_cambio_multiplicativo <- function(df) {
  bosque <- df %>% filter(pasture == 1)
  potrero <- df %>% filter(pasture == 0)
  
  abundancia_ratio <- bosque %>%
    select(starts_with("abun__draw_")) / potrero %>%
    select(starts_with("abun__draw_"))
  
  colnames(abundancia_ratio) <- gsub("abun__draw_", "ratio__draw", colnames(abundancia_ratio))
  
  df_resultado <- bosque %>%
    select(geometry, scientificName) %>%
    bind_cols(abundancia_ratio)
  
  return(df_resultado)
}

# function by ecoregion
sta_marta_resultado <- calcular_cambio_multiplicativo(sta_marta)
magdalena_dry_resultado <- calcular_cambio_multiplicativo(magdalena_dry)
cauca_montane_resultado <- calcular_cambio_multiplicativo(cauca_montane)
EC_montane_resultado <- calcular_cambio_multiplicativo(EC_montane)
villavicencio_dry_resultado <- calcular_cambio_multiplicativo(villavicencio_dry)

datos_ecoregiones <- bind_rows(
  mutate(sta_marta_resultado, ecorregion = "Santa Marta montane forests"),
  mutate(magdalena_dry_resultado, ecorregion = "Magdalena Valley dry forests"),
  mutate(cauca_montane_resultado, ecorregion = "Cauca Valley montane forests"),
  mutate(EC_montane_resultado, ecorregion = "Eastern Cordillera montane forests"),
  mutate(villavicencio_dry_resultado, ecorregion = "Villavicencio dry forests")
  )

library(tidyr)
datos_ecoregiones_limpios <- datos_ecoregiones %>%
  mutate(across(starts_with("ratio__draw"), ~ ifelse(is.finite(.), ., NA))) %>% 
  drop_na()

estadisticas1 <- datos_ecoregiones_limpios1 %>%
  pivot_longer(cols = starts_with("ratio__draw"),  # Convertir las columnas ratio__draw en formato largo
               names_to = "ratio", values_to = "valor") %>%
  group_by(ecorregion, ratio) %>%
  summarise(
    mediana = median(valor, na.rm = TRUE),
    p25 = quantile(valor, 0.25, na.rm = TRUE),
    p75 = quantile(valor, 0.75, na.rm = TRUE)
  )

library(ggridges)
ggplot(estadisticas1, aes(x = p25, y = reorder(ecorregion, p25), fill = ecorregion)) +
  geom_density_ridges(alpha = 0.6) +  
  labs(title = "25th percentile of the multiplicative abundance change across species",
       x = "sensitivity (N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none")

ggplot(estadisticas1, aes(x = mediana, y = reorder(ecorregion, mediana), fill = ecorregion)) +
  geom_density_ridges(alpha = 0.6) +  # Dibuja las áreas de densidad con un poco de transparencia
  labs(title = "Median multiplicative abundance change across species",
       x = "Sensitivity (N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none")

ggplot(estadisticas1, aes(x = p75, y = reorder(ecorregion, p75), fill = ecorregion)) +
  geom_density_ridges(alpha = 0.6) +  
  labs(title = "75th percentile of the multiplicative abundance change across species",
       x = "sensitivity (N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none")

ggplot(estadisticas1, aes(y = p25, x = reorder(ecorregion, -p25), fill = ecorregion)) +
  geom_boxplot(alpha = 0.6) +  
  labs(title = "25th percentile of the multiplicative abundance change across species",
       x = "Ecoregion",
       y = "Sensitivity (N forest/N pasture)") +
  theme_classic() +
  theme(legend.position = "none")

ggplot(estadisticas1, aes(y = mediana, x = reorder(ecorregion, -mediana), fill = ecorregion)) +
  geom_boxplot(alpha = 0.6) +  
  labs(title = "Median multiplicative abundance change across species",
       x = "Ecoregion",
       y = "Sensitivity (N forest/N pasture)") +
  theme_classic() +
  theme(legend.position = "none")

ggplot(estadisticas1, aes(y = p75, x = reorder(ecorregion, -p75), fill = ecorregion)) +
  geom_boxplot(alpha = 0.6) +  
  labs(title = "75th percentile of the multiplicative abundance change across species",
       x = "Ecoregion",
       y = "Sensitivity (N forest/N pasture)") +
  theme_classic() +
  theme(legend.position = "none")


col <- st_read("D:/Capas/Colombia/Colombia/COL_adm0.shp")
colombia_utm <- st_transform(col, crs = 32718)
grid_polygons <- st_make_grid(colombia_utm, cellsize = c(20000, 20000), what = "polygons")
grid_centroids <- st_centroid(grid_polygons)
plot(colombia_utm$geometry)
plot(grid_centroids, pch=20, col="blue", add=TRUE, cex=0.5)

###############################################################################
ecoregions_predictions <- readRDS("./Analysis/mean_abundance/ecoregions_predictions.rds")

multiplicative_change <- function(df) {
  forest <- df %>% filter(pasture == 1)
  pasture <- df %>% filter(pasture == 0)
  
  abundance_ratio <- forest %>%
    select(starts_with("abun__draw_")) / pasture %>%
    select(starts_with("abun__draw_"))
  
  colnames(abundance_ratio) <- gsub("abun__draw_", "ratio__draw_", colnames(abundance_ratio))
  
  result_df <- forest %>%
    select(scientificName) %>%
    bind_cols(abundance_ratio)
  return(result_df)
}

ratio_draw <- future_lapply(ecoregions_predictions, multiplicative_change, future.packages = c("dplyr", "sf"), future.seed = TRUE)
saveRDS(ratio_draw, "./ratio_draw_ecoregions.rds")

ratio_draw <- readRDS("./ratio_draw_ecoregions.rds")
mean_ratio <- function(df) {
  df %>%
    summarise(across(starts_with("ratio__draw_"), 
                     list(mean = ~mean(.x[is.finite(.x)], na.rm = TRUE),
                          p25  = ~quantile(.x[is.finite(.x)], probs = 0.25, na.rm = TRUE),
                          p50  = ~quantile(.x[is.finite(.x)], probs = 0.50, na.rm = TRUE),
                          p75  = ~quantile(.x[is.finite(.x)], probs = 0.75, na.rm = TRUE)))) %>% 
    mutate(ecoregion = df$ecorregion[1])
}

mean_ratio_draw <- bind_rows(lapply(ratio_draw, mean_ratio))

