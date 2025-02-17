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
  mutate(magdalena_dry_resultado, ecorregion = "Magdalena Valley dry forests")#,
#  mutate(cauca_montane_resultado, ecorregion = "Cauca Valley montane forests"),
#  mutate(EC_montane_resultado, ecorregion = "Eastern Cordillera montane forests"),
#  mutate(villavicencio_dry_resultado, ecorregion = "Villavicencio dry forests")
  )

library(tidyr)
datos_ecoregiones_limpios <- datos_ecoregiones %>%
  mutate(across(starts_with("ratio__draw"), ~ ifelse(is.finite(.), ., NA))) %>% 
  drop_na()

estadisticas <- datos_ecoregiones_limpios %>%
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

mean_ratio <- function(df, df_name) {
  df %>%
    summarise(across(starts_with("ratio__draw_"), 
                     list(mean = ~mean(.x[is.finite(.x)], na.rm = TRUE),
                          p25  = ~quantile(.x[is.finite(.x)], probs = 0.25, na.rm = TRUE),
                          p50  = ~quantile(.x[is.finite(.x)], probs = 0.50, na.rm = TRUE),
                          p75  = ~quantile(.x[is.finite(.x)], probs = 0.75, na.rm = TRUE)))) %>% 
    pivot_longer(cols = starts_with("ratio__draw_"),
                 names_to = c("ratio", ".value"),
                 names_pattern = "ratio__(draw_\\d+)_(.*)") %>%
    mutate(ecoregion = df_name)
}
mean_ratio_draw <- bind_rows(mapply(mean_ratio, ratio_draw, names(ratio_draw), SIMPLIFY = FALSE))

library(ggridges)
plot_p25 <- ggplot(mean_ratio_draw, aes(x = p25, y = reorder(ecoregion, p25), fill = ecoregion)) +
  geom_density_ridges2(alpha = 0.7) +  
  labs(title = "25th percentile of the multiplicative abundance change across species",
       x = "sensitivity (N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 8)) +
  stat_density(geom = "line", position = "identity", aes(x = p25), adjust = 1, kernel = "gaussian")

plot_mean <- ggplot(mean_ratio_draw, aes(x = mean, y = reorder(ecoregion, mean), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.6, scale = 1.5) +  
  labs(title = "Mean of the multiplicative abundance change across species",
       x = "sensitivity (N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 8)) +
  stat_density(geom = "line", position = "identity", aes(x = mean), adjust = 1, kernel = "gaussian")

plot_p50 <- ggplot(mean_ratio_draw, aes(x = p50, y = reorder(ecoregion, p50), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.6, scale = 1.5) +  
  labs(title = "Median of the multiplicative abundance change across species",
       x = "sensitivity (N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 8)) +
  stat_density(geom = "line", position = "identity", aes(x = p50), adjust = 1, kernel = "gaussian")

plot_p75 <- ggplot(mean_ratio_draw, aes(x = p75, y = reorder(ecoregion, p75), fill = ecoregion)) +
  geom_density_ridges2(alpha = 0.6) +  
  labs(title = "75th percentile of the multiplicative abundance change across species",
       x = "sensitivity (N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 8)) +
  stat_density(geom = "line", position = "identity", aes(x = p75), kernel = "gaussian")

library(gridExtra)

grid.arrange(plot_p25, plot_p50, plot_p75,  ncol=1)

#############################################
ecoregions_predictions <- readRDS("./Analysis/mean_abundance/ecoregions_predictions.rds")
#draft <- ecoregions_predictions[c(8, 11, 13)]
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(sf)

# Get the region names from the list names
region_names <- gsub(" forests", "", names(ecoregions_predictions))

# Function to calculate the proportional change in abundance by species and region
abundance_change <- function(sf_df, region_name) {
  sf_df %>%
    st_drop_geometry() %>% 
    # Divide into forest (pasture = 1) and grassland (pasture = 0)
    group_by(scientificName, pasture) %>%
    summarise(across(starts_with("abun__draw_"), \(x) sum(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_wider(names_from = pasture, values_from = starts_with("abun__draw_"), 
                names_glue = "{.value}_pasture{pasture}") %>%
    mutate(across(ends_with("_pasture1"), 
                  ~ .x / get(sub("_pasture1", "_pasture0", cur_column())), 
                  .names = "ratio_{.col}")) %>%
    select(scientificName, starts_with("ratio_")) %>%
    mutate(region = region_name)
}

# Apply the function to each element of the list
ratio_draw <- map2(ecoregions_predictions, region_names, abundance_change)

ratio_names <- gsub(" forests", "", names(ratio_draw))

mean_ratio <- function(df, df_name) {
  df %>%
    summarise(across(starts_with("ratio_abun__draw_"), 
                     list(mean = ~mean(.x[is.finite(.x)], na.rm = TRUE),
                          p25  = ~quantile(.x[is.finite(.x)], probs = 0.25, na.rm = TRUE),
                          p50  = ~quantile(.x[is.finite(.x)], probs = 0.50, na.rm = TRUE),
                          p75  = ~quantile(.x[is.finite(.x)], probs = 0.75, na.rm = TRUE)))) %>% 
    pivot_longer(cols = starts_with("ratio_abun__draw_"),
                 names_to = c("draw", ".value"),
                 names_pattern = "ratio_abun__draw_(\\d+)_(.*)") %>%
    mutate(ecoregion = df_name)
}
mean_ratio_draw <- bind_rows(mapply(mean_ratio, ratio_draw, names(ratio_draw), SIMPLIFY = FALSE))
mean_ratio_draw$ecoregion <- gsub(" forests", "", mean_ratio_draw$ecoregion)

library(ggridges)
plot_p25 <- ggplot(mean_ratio_draw, aes(x = pasture1_p25, y = reorder(ecoregion, -pasture1_p25), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.7) +  
  labs(title = "25th percentile",
       x = "sensitivity (N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=6)) +
  scale_x_continuous(limits = c(0, 8)) +
  stat_density(geom = "line", position = "identity", aes(x = pasture1_p25), adjust = 1, kernel = "gaussian")

plot_mean <- ggplot(mean_ratio_draw, aes(x = pasture1_mean, y = reorder(ecoregion, -pasture1_mean), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.6, scale = 1.5) +  
  labs(title = "Mean of the multiplicative abundance change across species",
       x = "sensitivity (N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 8)) +
  stat_density(geom = "line", position = "identity", aes(x = pasture1_mean), adjust = 1, kernel = "gaussian")

plot_p50 <- ggplot(mean_ratio_draw, aes(x = pasture1_p50, y = reorder(ecoregion, -pasture1_p50), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.6, scale = 1.5) +  
  labs(title = "Median",
       x = "sensitivity (N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=6)) +
  scale_x_continuous(limits = c(0, 8)) +
  stat_density(geom = "line", position = "identity", aes(x = pasture1_p50), adjust = 1, kernel = "gaussian")

plot_p75 <- ggplot(mean_ratio_draw, aes(x = pasture1_p75, y = reorder(ecoregion, -pasture1_p75), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.6) +  
  labs(title = "75th percentile",
       x = "sensitivity (N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=6)) +
  scale_x_continuous(limits = c(0, 8)) +
  stat_density(geom = "line", position = "identity", aes(x = pasture1_p75), kernel = "gaussian")

#library(gridExtra)

grid.arrange(plot_p25, plot_p50, plot_p75,  ncol=1)

marta <- ecoregions_predictions[[13]]
marta_forest <- subset(marta, marta$pasture==1)  
marta_pasture <- subset(marta, marta$pasture==0)  

sum(marta_forest$abun__draw_300[marta_forest$scientificName == "Canthidium_sp._01H"], na.rm = TRUE) /
  sum(marta_pasture$abun__draw_300[marta_pasture$scientificName == "Canthidium_sp._01H"], na.rm = TRUE)

