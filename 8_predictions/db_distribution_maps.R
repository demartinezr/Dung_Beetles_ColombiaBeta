setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata")

library(sf)
library(dplyr)
library(ggplot2)

# distribution maps for dung beetles abudance across Colombia
#
# abundance predictions for 243 species of dung beetles with 10 iterations in 
# pasture and forest
  db_predictions <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/species_predictions.rds")
# study area ecoregions and grid
  study_area <- st_read("D:/Capas/America/ecoregions/ecoreg.shp")
  grid_2km <- st_read("D:/Capas/America/grid/grid_2km_ID_WGS84.shp")
#
# mean abundance for Ateuchus cracicus an endemic species from Magdalena Valley
  
  cracicus <- db_predictions[["Ateuchus_cracicus"]]
    cracicus <- cracicus %>% select(-abun__draw_600)
    cracicus_pasture <- subset(cracicus, cracicus$pasture == 0)
    cracicus_forest <- subset(cracicus, cracicus$pasture == 1)

  # Calculate mean abundance in forest and pasture scenarios across Colombia
  forest_mean <- cracicus_forest %>%
    rowwise() %>%
    mutate(mean_abundance = mean(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
    ungroup()

  pasture_mean <- cracicus_pasture %>% 
    rowwise() %>%
    mutate(mean_abundance = mean(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
    ungroup()
  # join abundance mean to grid
    grid_forest <- grid_2km %>%
    st_join(forest_mean, left = TRUE)
  
  grid_forest <- grid_forest %>%
    mutate(mean_abundance = ifelse(is.na(mean_abundance), 0, mean_abundance))
  # distribution map in forest and pasture
  ggplot(grid_forest) +
    geom_sf(aes(fill = mean_abundance), color = NA) + 
    scale_fill_gradient(low = "blue", high = "darkblue", name = "Abundance Mean", labels = scales::label_number()) +
    theme_classic() +
    labs(title = "Ateuchus cracicus forest")
  
  grid_pasture <- grid_2km %>%
    st_join(pasture_mean, left = TRUE)
  
  grid_pasture <- grid_pasture %>%
    mutate(mean_abundance = ifelse(is.na(mean_abundance), 0, mean_abundance))
  
  ggplot(grid_pasture) +
    geom_sf(aes(fill = mean_abundance), color = NA) + 
    scale_fill_gradient(low = "blue", high = "darkblue", name = "Abundance Mean") +
    theme_classic() +
    labs(title = "Ateuchus cracicus pasture")
#
  # mean abundance for Digitonthophagus gazella an introduced species
  
  gazella <- db_predictions[["Digitonthophagus_gazella"]]
  gazella <- gazella %>% select(-abun__draw_600)
  gazella_pasture <- subset(gazella, gazella$pasture == 0)
  gazella_forest <- subset(gazella, gazella$pasture == 1)
  
  # Calculate median abundance in forest and pasture scenarios across Colombia
  forest_median <- gazella_forest %>%
    rowwise() %>%
    mutate(median_abundance = median(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
    ungroup()
  
  pasture_median <- gazella_pasture %>% 
    rowwise() %>%
    mutate(median_abundance = median(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
    ungroup()
  # join abundance mean to grid
  grid_forest <- grid_2km %>%
    st_join(forest_median, left = TRUE)
  
  grid_forest <- grid_forest %>%
    mutate(median_abundance = ifelse(is.na(median_abundance), 0, median_abundance))
  # distribution map in forest and pasture
  ggplot(grid_forest) +
    geom_sf(aes(fill = median_abundance), color = NA) + 
    scale_fill_gradient(low = "blue", high = "red", name = "Median abundance", labels = scales::label_number()) +
    theme_classic() +
    labs(title = "Digitonthophagus gazella forest")
  
  grid_pasture <- grid_2km %>%
    st_join(pasture_median, left = TRUE)
  
  grid_pasture <- grid_pasture %>%
    mutate(median_abundance = ifelse(is.na(median_abundance), 0, median_abundance))
  
  ggplot(grid_pasture) +
    geom_sf(aes(fill = median_abundance), color = NA) + 
    scale_fill_viridis_c(name = "Median abundance", option = "B") +
    theme_classic() +
    labs(title = "Digitonthophagus gazella pasture")
  
  sum(gazella$pasture == 0)
  
  ggplot(data = gazella_pasture) +
    geom_sf(aes(color = abun__draw_1200), size = 3) +
    scale_color_viridis_c(name = "Median Abundance", option = "D") +
    theme_minimal() +
    labs(title = "D. gazella median abundance in pasture")
  
  # # mean abundance for Digitonthophagus gazella an introduced species
  
  marmoreus <- db_predictions[["Eurysternus_marmoreus"]]
#  gazella <- gazella %>% select(-abun__draw_600)
  marmoreus_pasture <- subset(marmoreus, marmoreus$pasture == 0)
  marmoreus_forest <- subset(marmoreus, marmoreus$pasture == 1)
  
  # Calculate median abundance in forest and pasture scenarios across Colombia
  forest_median <- marmoreus_forest %>%
    rowwise() %>%
    mutate(median_abundance = median(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
    ungroup()
  
  pasture_median <- marmoreus_pasture %>% 
    rowwise() %>%
    mutate(median_abundance = median(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
    ungroup()
  # join abundance mean to grid
  grid_forest <- grid_2km %>%
    st_join(forest_median, left = TRUE)
  
  grid_forest <- grid_forest %>%
    mutate(median_abundance = ifelse(is.na(median_abundance), 0, median_abundance))
  # distribution map in forest and pasture
  ggplot(grid_forest) +
    geom_sf(aes(fill = median_abundance), color = NA) + 
    scale_color_viridis_c(name = "Median Abundance", option = "D") +
    theme_classic() +
    labs(title = "Eurysternus marmoreus forest")
  
  grid_pasture <- grid_2km %>%
    st_join(pasture_median, left = TRUE)
  
  grid_pasture <- grid_pasture %>%
    mutate(median_abundance = ifelse(is.na(median_abundance), 0, median_abundance))
  
  ggplot(grid_pasture) +
    geom_sf(aes(fill = median_abundance), color = NA) + 
    scale_color_viridis_c(name = "Median Abundance", option = "D") +
    theme_classic() +
    labs(title = "Eurysternus marmoreus pasture")
  
  sum(gazella$pasture == 0)
  
  ggplot(data = marmoreus_pasture) +
    geom_sf(aes(color = abun__draw_1200), size = 3) +
    scale_color_viridis_c(name = "Median Abundance", option = "D") +
    theme_minimal() +
    labs(title = "E. marmoreus median abundance in pasture")
  