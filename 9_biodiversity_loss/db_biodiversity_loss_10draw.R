setwd("C:/Users/PC/Dropbox/CO_DBdata")
ecoregions_predictions <- readRDS("./Analysis/mean_abundance/ecoregions_predictions_2.rds")

library(raster)
library(purrr)
library(dplyr)
library(data.table)
library(ggplot2)

AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# ecoregions names in each sf dataframe
ecoregions_predictions <- map2(
  ecoregions_predictions, 
  names(ecoregions_predictions), 
  ~ mutate(.x, ecoregions = .y)
)
# from list to sf dataframe
ecoregions_all <- reduce(ecoregions_predictions, rbind)
#saveRDS(ecoregions_all, "./ecoregions_all_2.rds")
# filter draws
ecoregions_all <- readRDS("./ecoregions_all.rds")

# load functions to compute biodiversity loss at local and reigonal scales for all draws
source("F:/repositorio/Dung_Beetles_ColombiaBeta/9_biodiversity_loss/db_compute_loss.R")

# Draws
#draws <- seq(300, 3000, by = 300)
draws <-  setdiff(seq(100, 3000, by = 100), 30 * (1:100))[1:20]

# List
cell_ratios_list <- list()
colombia_ratios_list <- list()

gc()  

for (i in seq_along(draws)) {
  cat("draw:", i, "of", length(draws), "\n")  
  draw_col <- paste0("abun__draw_", draws[i])
  
  # Filter data and add ID cell
  data_draw <- ecoregions_all %>%
    select(ecoregions, lon, lat, scientificName, pasture, all_of(draw_col)) %>%
    mutate(cell_id = row_number()) %>%
    select(cell_id, everything())
  row.names(data_draw) <- NULL
  
  # forest and pasture data sets
  forest_draw <- data_draw %>% filter(pasture == 1) %>% select(-pasture)
  pasture_draw <- data_draw %>% filter(pasture == 0) %>% select(-pasture)
  
  # Convertir a formato largo y transformar en tablas anchas
  forest_draw_dt <- as.data.table(forest_draw)
  forest_draw_wide <- dcast(forest_draw_dt, ecoregions + lon + lat ~ scientificName, 
                            value.var = draw_col, fill = 0)
  
  pasture_draw_dt <- as.data.table(pasture_draw)
  pasture_draw_wide <- dcast(pasture_draw_dt, ecoregions + lon + lat ~ scientificName, 
                             value.var = draw_col, fill = 0)
  
  # Calcular métricas
  cell_ratios_list[[draw_col]] <- get_avg_cell_ratios(forest_draw_wide, pasture_draw_wide, 
                                                      cutoff_type = "relative", cutoff = .2)
  colombia_ratios_list[[draw_col]] <- get_regional_ratios(forest_draw_wide, pasture_draw_wide, 
                                                          cutoff_type = "relative", cutoff = .2, 
                                                          cell_positions = NULL)
  
  # memory check
  rm(forest_draw, pasture_draw, forest_draw_dt, forest_draw_wide, pasture_draw_dt, pasture_draw_wide, data_draw)
  gc()
}

# Combinar resultados
cell_ratios_df <- bind_rows(cell_ratios_list, .id = "draw")
colombia_ratios_df <- bind_rows(colombia_ratios_list, .id = "draw")

saveRDS(cell_ratios_df, "./cell_ratios_df_2.rds")
saveRDS(colombia_ratios_df, "./colombia_ratios_df_2.rds")


##### Colombia-wide comparison #####
mean(cell_ratios_df$med_logratio, na.rm = T)

mean(exp(cell_ratios_df$med_logratio), na.rm = T)
exp(mean(cell_ratios_df$med_logratio, na.rm = T))
exp(colombia_ratios_df$med_logratio)

# Calcular promedios de métricas
mean(cell_ratios_df$med_logratio, na.rm = TRUE)
mean(exp(cell_ratios_df$med_logratio), na.rm = TRUE)
exp(mean(cell_ratios_df$med_logratio, na.rm = TRUE))

mean(unlist(lapply(colombia_ratios_list, function(x) x$med_logratio)), na.rm = TRUE)
exp(mean(unlist(lapply(colombia_ratios_list, function(x) x$med_logratio)), na.rm = TRUE))

# Mapeo de log-ratios
ggplot(cell_ratios_df, aes(x = lon, y = lat, fill = avg_logratio)) +
  geom_raster() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(title = "Local scale diversity loss", fill = "Log ratio")

# Mapeo de riqueza de especies
ggplot(cell_ratios_df, aes(x = lon, y = lat, fill = n)) +
  geom_raster() +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(title = "Species richness", fill = "Number of species")

#################### ecoregions analysis #####################################

library(ggplot2)
library(ggridges)
library(gridExtra)

ecoregions_all <- readRDS("./ecoregions_all.rds")

source("F:/repositorio/Dung_Beetles_ColombiaBeta/9_biodiversity_loss/db_compute_loss.R")
draws <- seq(300, 3000, by = 300)
#draws <-  setdiff(seq(100, 3000, by = 100), 30 * (1:100))[1:20]
unique_ecoregions <- unique(ecoregions_all$ecoregions)

cell_ratios_list <- list()
colombia_ratios_list <- list()
ecoregion_ratios_list <- list()

for (i in seq_along(draws)) {
  cat("draw:", i, "of", length(draws), "\n")  
  draw_col <- paste0("abun__draw_", draws[i])
  
  # Filter data and add ID cell
  data_draw <- ecoregions_all %>%
    select(ecoregions, lon, lat, scientificName, pasture, all_of(draw_col)) %>%
    mutate(cell_id = row_number()) %>%
    select(cell_id, everything())
  row.names(data_draw) <- NULL
  
  # forest and pasture data sets
  forest_draw <- data_draw %>% filter(pasture == 0) %>% select(-pasture)
  pasture_draw <- data_draw %>% filter(pasture == 1) %>% select(-pasture)
  
  # format data sets
  forest_draw_dt <- as.data.table(forest_draw)
  forest_draw_wide <- dcast(forest_draw_dt, ecoregions + lon + lat ~ scientificName, 
                            value.var = draw_col, fill = 0)
  
  pasture_draw_dt <- as.data.table(pasture_draw)
  pasture_draw_wide <- dcast(pasture_draw_dt, ecoregions + lon + lat ~ scientificName, 
                             value.var = draw_col, fill = 0)
  
  cell_ratios_list[[draw_col]] <- get_avg_cell_ratios(forest_draw_wide, pasture_draw_wide, 
                                                      cutoff_type = "absolute", cutoff = 1)
  colombia_ratios_list[[draw_col]] <- get_regional_ratios(forest_draw_wide, pasture_draw_wide, 
                                                          cutoff_type = "absolute", cutoff = 1, 
                                                          cell_positions = NULL)

  ecoregion_results <- list()
  
  for (j in seq_along(unique_ecoregions)) {
    cat(j, "\n") 
    
    eco <- unique_ecoregions[j]
    cell_positions <- which(forest_draw_wide$ecoregions == eco)
    eco_result <- get_regional_ratios(forest_draw_wide, pasture_draw_wide, 
                                      cutoff_type = "absolute", cutoff = 1, 
                                      cell_positions = cell_positions)
    
    eco_result$ecoregion <- eco  # Agregar la ecorregión a los resultados
    ecoregion_results[[eco]] <- eco_result
  }
  
  ecoregion_ratios_list[[draw_col]] <- ecoregion_results
  rm(forest_draw, forest_draw_dt, forest_draw_wide, 
     pasture_draw, pasture_draw_dt, pasture_draw_wide, data_draw)
  gc()
}

saveRDS(cell_ratios_list, "./diversity_loss/cell_ratios_list.rds")
saveRDS(ecoregion_ratios_list, "./diversity_loss/ecoregion_ratios_list.rds")
saveRDS(colombia_ratios_list, "./diversity_loss/colombia_ratios_list.rds")

######################### ecoregions dataframe #################################
#cell_ratios_r01 <- readRDS("./cell_ratios_list_r01.rds")
#cell_ratios_r01df <- bind_rows(cell_ratios_r01, .id = "draw_col")

#ecoregion_ratios_c02 <- readRDS("./ecoregion_ratios_list_c02.rds")
#ecoregion_ratios_02df <- bind_rows(lapply(ecoregion_ratios_c02, function(draw) {
#  bind_rows(draw, .id = "ecoregion")
#}), .id = "draw_col")

#colombia_ratios_20draw <- readRDS("./colombia_ratios_list_0.2_20draw.rds")
#colombia_ratios_20drawdf <- bind_rows(colombia_ratios_20draw, .id = "draw_col")
#colombia_ratios_20drawdf$ecoregion <- "Colombia"


cell_ratios_list <- readRDS("./diversity_loss/cell_ratios_list.rds")
cell_ratios_df <- bind_rows(cell_ratios_list, .id = "draw_col")

ecoregion_ratios_list <- readRDS("./diversity_loss/ecoregion_ratios_list.rds")
ecoregion_ratios_df <- bind_rows(lapply(ecoregion_ratios_list, function(draw) {
  bind_rows(draw, .id = "ecoregion")
}), .id = "draw_col")

colombia_ratios_list <- readRDS("./diversity_loss/colombia_ratios_list_c02.rds")
colombia_ratios_df <- bind_rows(lapply(colombia_ratios_list, function(draw) {
  bind_rows(draw) %>% mutate(ecoregion = "Colombia")
}), .id = "draw_col")


##### Colombia-wide comparison 10 draws #####
mean(cell_ratios_df$med_logratio, na.rm = T)
mean(exp(cell_ratios_df$med_logratio), na.rm = T)
exp(mean(cell_ratios_df$med_logratio, na.rm = T))

mean(colombia_ratios_df$med_logratio, na.rm = T)
mean(exp(colombia_ratios_df$med_logratio))
exp(mean(colombia_ratios_df$med_logratio))


##### Colombia-wide comparison 20 draws #####
mean(cell_ratios_df$med_logratio, na.rm = T)
mean(exp(cell_ratios_df$med_logratio), na.rm = T)
exp(mean(cell_ratios_df$med_logratio, na.rm = T))

mean(colombia_ratios_df$med_logratio, na.rm = T)
mean(exp(colombia_ratios_df$med_logratio))
exp(mean(colombia_ratios_df$med_logratio))

plot_df <- bind_rows(ecoregion_ratios_df, colombia_ratios_df)
plot_df$ecoregion <- gsub(" forests", "", plot_df$ecoregion)

# plot_20drawdf <- bind_rows(ecoregion_ratios_20drawdf, colombia_ratios_20drawdf)
# plot_20drawdf$ecoregion <- gsub(" forests", "", plot_20drawdf$ecoregion)

p1 <- ggplot(plot_df, aes(x = reorder(ecoregion, -p_25_logratio), y = p_25_logratio, fill = ecoregion)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(title = "25th percentile",
       x = NULL,
       y = "Logratio") +
  theme_bw() +
#  scale_y_continuous(limits=c(-8,1)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8))

p2 <- ggplot(plot_df, aes(x = reorder(ecoregion, -avg_logratio), y = avg_logratio, fill = ecoregion)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(title = "Mean",
       x = NULL,
       y = "Logratio") +
  theme_bw() +
#  scale_y_continuous(limits=c(-8,1)) +
  theme(legend.position = "none",axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8))

p3 <- ggplot(plot_df, aes(x = reorder(ecoregion, -p_75_logratio), y = p_75_logratio, fill = ecoregion)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(title = "75th percentile",
       x = NULL,
       y = "Logratio") +
  theme_bw() +
#  scale_y_continuous(limits=c(-8,1)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8))


grid.arrange(p1, p2, p3, ncol=3)

p1_10 <- ggplot(plot_20drawdf, aes(x = reorder(ecoregion, -p_25_logratio), y = p_25_logratio, fill = ecoregion)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(title = "25th percentile",
       x = NULL,
       y = "Logratio") +
  theme_bw() +
#  scale_y_continuous(limits=c(-3,1)) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8))

p2_10 <- ggplot(plot_20drawdf, aes(x = reorder(ecoregion, -avg_logratio), y = avg_logratio, fill = ecoregion)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(title = "Mean",
       x = NULL,
       y = "Logratio") +
  theme_bw() +
#  scale_y_continuous(limits=c(-3,1)) +
  theme(legend.position = "none",axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8))

p3_10 <- ggplot(plot_20drawdf, aes(x = reorder(ecoregion, -p_75_logratio), y = p_75_logratio, fill = ecoregion)) +
  geom_boxplot(color = "black", alpha = 0.7) +
  scale_fill_viridis_d() +
  labs(title = "75th percentile",
       x = NULL,
       y = "Logratio") +
  theme_bw() +
#  scale_y_continuous(limits=c(-3,1)) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8))

grid.arrange(p1_10, p2_10, p3_10, p1, p2, p3, ncol=3)

################## ecoregions comparison #################
plot_mean_df <- plot_df %>%
  group_by(ecoregion) %>%
  summarise(
    avg_ratio = mean(avg_ratio, na.rm = TRUE),
    avg_logratio = mean(avg_logratio, na.rm = TRUE),
    med_logratio = mean(med_logratio, na.rm = TRUE),
    p_25_logratio = mean(p_25_logratio, na.rm = TRUE),
    p_75_logratio = mean(p_75_logratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(ecoregion != "Colombia") %>% 
  arrange(avg_logratio) %>%
  mutate(ecoregion_order = row_number())

######## ecoregions diversity loss across the posterior distribution ###########


plot_mean <- ggplot(plot_df, aes(x = avg_logratio, y = reorder(ecoregion, -avg_logratio), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5)  +
  scale_fill_viridis_d() +
  labs(title = "Mean", x = "sensitivity\n(N forest/N pasture)", y = "Ecoregion") +
  theme_bw() +
  theme(legend.position = "none",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits=c(-3, 4)) +
  stat_density(geom = "line", position = "identity", aes(x = avg_logratio), kernel = "gaussian")


plot_25 <- ggplot(plot_df, aes(x = p_25_logratio, y = reorder(ecoregion, -p_25_logratio), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5)  +
  scale_fill_viridis_d() +
  labs(title = "25th percentile", x = "sensitivity\n(N forest/N pasture)", y = "Ecoregion") +
  theme_bw() +
  theme(legend.position = "none",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 10, face = "bold")) +
    scale_x_continuous(limits=c(-3, 4)) +
  stat_density(geom = "line",position = "identity", aes(x = p_25_logratio), kernel = "gaussian")

plot_75 <- ggplot(plot_df, aes(x = p_75_logratio, y = reorder(ecoregion, -p_75_logratio), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5)  +
  scale_fill_viridis_d() +
  labs(title = "75th percentile", x = "sensitivity (N forest/N pasture)", y = "Ecoregion") +
  theme_bw() +
  theme(legend.position = "none",
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits=c(-3, 4)) +
  stat_density(geom = "line", position = "identity", aes(x = p_75_logratio), kernel = "gaussian")

grid.arrange(plot_25, plot_mean, plot_75, ncol = 1)

###################### Relative difference #############################
relative_diff <- plot_df %>%
  left_join(colombia_ratios_df, by = "draw_col", suffix = c("", "_col")) %>% 
  mutate(
    avg_ratio_dif = avg_ratio / avg_ratio_col,
    avg_logratio_dif  = exp(avg_logratio_col - avg_logratio),
    med_logratio_dif  = exp(med_logratio_col - med_logratio),
    p_25_logratio_dif  = exp(p_25_logratio_col - p_25_logratio),
    p_75_logratio_dif  = exp(p_75_logratio_col - p_75_logratio)) %>%
  select(ecoregion, draw_col, avg_ratio_dif, avg_logratio_dif, med_logratio_dif, p_25_logratio_dif, p_75_logratio_dif)


relative_diff_mean <-relative_diff %>%
  group_by(ecoregion) %>%
  summarise(
    avg_ratio_dif = mean(avg_ratio_dif, na.rm = TRUE),
    avg_logratio_dif = mean(avg_logratio_dif, na.rm = TRUE),
    med_logratio_dif = mean(med_logratio_dif, na.rm = TRUE),
    p_25_logratio_dif = mean(p_25_logratio_dif, na.rm = TRUE),
    p_75_logratio_dif = mean(p_75_logratio_dif, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(ecoregion != "Colombia") %>% 
  arrange(avg_logratio_dif) %>%
  mutate(ecoregion_order = row_number())


plot_mean_dif <- ggplot(relative_diff, aes(x = avg_logratio_dif, y = reorder(ecoregion, avg_logratio_dif), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5)  +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
  scale_fill_viridis_d() +
  labs(title = "Mean", x = "sensitivity\n(N forest/N pasture)", y = "Ecoregion") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits = c(0, 5)) +
  stat_density(geom = "line", position = "identity", aes(x = avg_logratio_dif), kernel = "gaussian")


plot_25_dif <- ggplot(relative_diff, aes(x = p_25_logratio_dif, y = reorder(ecoregion, p_25_logratio_dif), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5)  +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
  scale_fill_viridis_d() +
  labs(title = "25th percentile", x = "sensitivity\n(N forest/N pasture)", y = "Ecoregion") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits = c(0, 5)) +
  stat_density(geom = "line",position = "identity", aes(x = p_25_logratio_dif), kernel = "gaussian")

plot_75_dif <- ggplot(relative_diff, aes(x = p_75_logratio_dif, y = reorder(ecoregion, p_75_logratio_dif), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5)  +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
  scale_fill_viridis_d() +
  labs(title = "75th percentile", x = "Relative difference", y = "Ecoregion") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits = c(0, 5)) +
  stat_density(geom = "line", position = "identity", aes(x = p_75_logratio_dif), kernel = "gaussian")

grid.arrange(plot_25, plot_25_dif, plot_mean, plot_mean_dif, plot_75, plot_75_dif, nrow = 3, widths =c(1.5,1))


################ Relative difference cummulative ####################################
library(tidyr)
plot_mean_df <- plot_df %>%
  group_by(ecoregion) %>%
  summarise(
    avg_ratio = mean(avg_ratio, na.rm = TRUE),
    avg_logratio = mean(avg_logratio, na.rm = TRUE),
    med_logratio = mean(med_logratio, na.rm = TRUE),
    p_25_logratio = mean(p_25_logratio, na.rm = TRUE),
    p_75_logratio = mean(p_75_logratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(ecoregion != "Colombia") %>% 
  arrange(avg_logratio) %>%
  mutate(ecoregion_order = row_number())

pan_colombia_avg <- colombia_ratios_df %>%
  summarise(
    avg_logratio_mean  = mean(avg_logratio, na.rm = TRUE),
    med_logratio_mean  = mean(med_logratio, na.rm = TRUE),
    p25_logratio_mean  = mean(p_25_logratio, na.rm = TRUE),
    p75_logratio_mean  = mean(p_75_logratio, na.rm = TRUE),
    .groups = "drop"
  )

plot_mean_df <- plot_mean_df %>%
  mutate(relative_difference = pan_colombia_avg$med_logratio_mean - med_logratio)

set.seed(123) # Para reproducibilidad

num_sequences <- 1000
num_ecoregions <- nrow(plot_mean_df) # Cantidad de ecorregiones en el análisis

generate_sensitivity_sequence <- function() {
  shuffled_ecoregions <- sample(plot_mean_df$ecoregion)  # Mezclar aleatoriamente las ecorregiones
  
  temp_df <- plot_mean_df %>%
    filter(ecoregion %in% shuffled_ecoregions) %>%
    arrange(match(ecoregion, shuffled_ecoregions)) %>%
    mutate(
      num_regions = row_number(),
      cum_avg_logratio  = cumsum(avg_logratio),
      cum_med_logratio  = cumsum(med_logratio),
      cum_p25_logratio  = cumsum(p_25_logratio),
      cum_p75_logratio  = cumsum(p_75_logratio),
      relative_diff_avg = pan_colombia_avg$avg_logratio_mean - cum_avg_logratio,
      relative_diff_med = pan_colombia_avg$med_logratio_mean - cum_med_logratio,
      relative_diff_p25 = pan_colombia_avg$p25_logratio_mean - cum_p25_logratio,
      relative_diff_p75 = pan_colombia_avg$p75_logratio_mean - cum_p75_logratio
    ) %>%
    select(num_regions, relative_diff_avg, relative_diff_med, relative_diff_p25, relative_diff_p75) %>%
    drop_na()  # Eliminar valores problemáticos
  
  return(temp_df)
}

# Generar 1000 simulaciones
sensitivity_results <- replicate(1000, generate_sensitivity_sequence(), simplify = FALSE)

# Convertir a un solo data frame
sensitivity_df <- bind_rows(sensitivity_results, .id = "simulation_id")

# Crear `sensitivity_summary`
sensitivity_summary <- sensitivity_df %>%
  group_by(num_regions) %>%
  summarise(
    mean_relative_diff_avg = mean(relative_diff_avg),
    lower_ci_avg = quantile(relative_diff_avg, probs = 0.1),
    upper_ci_avg = quantile(relative_diff_avg, probs = 0.9),
    mean_relative_diff_med = mean(relative_diff_med),
    lower_ci_med = quantile(relative_diff_med, probs = 0.1),
    upper_ci_med = quantile(relative_diff_med, probs = 0.9),
    mean_relative_diff_p25 = mean(relative_diff_p25),
    lower_ci_p25 = quantile(relative_diff_p25, probs = 0.1),
    upper_ci_p25 = quantile(relative_diff_p25, probs = 0.9),
    mean_relative_diff_p75 = mean(relative_diff_p75),
    lower_ci_p75 = quantile(relative_diff_p75, probs = 0.1),
    upper_ci_p75 = quantile(relative_diff_p75, probs = 0.9),
    .groups = "drop"
  )


# Visualización del resultado
library(ggplot2)

plot_mean_cum <- ggplot(sensitivity_summary, aes(x = num_regions, y = exp(mean_relative_diff_avg))) +
  geom_point(size = 3, color = "black") +  # Puntos para la media
  geom_errorbar(aes(ymin = exp(lower_ci_avg), ymax = exp(upper_ci_avg)), width = 0.2, color = "black") +  # Barras de error
  geom_ribbon(aes(ymin = -0.05, ymax = 0.05), fill = "gray", alpha = 0.5) +  # Banda gris
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Línea de referencia en 1
  scale_x_continuous(breaks = seq(1, 13, by = 2)) +  # Escala del eje X
#  scale_y_continuous(limits=c(0,8)) +
  labs(title = "Mean",
    x = NULL,
    y = "Relative difference"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10))# +
 # scale_y_continuous(limits = c(0, 6)))

plot_25_cum <- ggplot(sensitivity_summary, aes(x = num_regions, y = exp(mean_relative_diff_p25))) +
  geom_point(size = 3, color = "black") +  # Puntos para la media
  geom_errorbar(aes(ymin = exp(lower_ci_p25), ymax = exp(upper_ci_p25)), width = 0.2, color = "black") +  # Barras de error
  geom_ribbon(aes(ymin = -0.05, ymax = 0.05), fill = "gray", alpha = 0.5) +  # Banda gris
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Línea de referencia en 1
  scale_x_continuous(breaks = seq(1, 13, by = 2)) +  # Escala del eje X
#  scale_y_continuous(limits=c(0,8)) +
  labs(title="25th percentile",
    x = NULL,
    y = "Relative difference"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10)
  )

plot_75_cum <- ggplot(sensitivity_summary, aes(x = num_regions, y = exp(mean_relative_diff_p75))) +
  geom_point(size = 3, color = "black") +  # Puntos para la media
  geom_errorbar(aes(ymin = exp(lower_ci_p75), ymax = exp(upper_ci_p75)), width = 0.2, color = "black") +  # Barras de error
  geom_ribbon(aes(ymin = -0.05, ymax = 0.05), fill = "gray", alpha = 0.5) +  # Banda gris
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Línea de referencia en 1
  scale_x_continuous(breaks = seq(1, 13, by = 2)) +  # Escala del eje X
  scale_y_continuous(limits=c(0,8)) +
  labs(title ="75th percentile",
    x = expression(N[regions]),
    y = "Relative difference"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size  = 14, face = "italic"),
    axis.title.y = element_text(size  = 10)
  )

grid.arrange(plot_25, plot_25_dif, plot_25_cum, 
             plot_mean, plot_mean_dif, plot_mean_cum,
             plot_75, plot_75_dif, plot_75_cum, ncol= 3, widths =c(1.5,1,1))

grid.arrange(plot_25_cum, plot_mean_cum, plot_75_cum)

################################################################################
####################### Grids with varying pixel size ##########################
library(dggridR)
library(sf)
library(raster)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggridges)
library(gridExtra)
setwd("C:/Users/PC/Dropbox/CO_DBdata")

ecoregions_all <- as.data.frame(readRDS("./ecoregions_all.rds"))

source("F:/repositorio/Dung_Beetles_ColombiaBeta/9_biodiversity_loss/db_compute_loss.R")


# Assign each point to a hexagonal cell on a grid of variable resolution
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
forest_data <- ecoregions_all %>% filter(pasture==1)
points_latlong <- forest_data %>%  select(ecoregions, lat, lon) %>%  distinct() %>%  st_as_sf(coords = c("lon", "lat"), crs = 4326)
points_coords <- as.data.frame(st_coordinates(points_latlong))

# Definir resoluciones y draws
resolutions <- c(5:8)
draws <- seq(300, 3000, by = 300)
draws <- c(300, 1500, 3000)

# Inicializar listas para almacenar datos
forest_list <- vector("list", length(draws))
pasture_list <- vector("list", length(draws))
names(forest_list) <- names(pasture_list) <- paste0("draw_", draws)

# Generar forest_*_wide y pasture_*_wide para cada draw
for (draw in draws) {
  cat("draw", draw, "\n")
  draw_col <- paste0("abun__draw_", draw)
  
  ecoregions_tmp <- ecoregions_all %>%
    select(ecoregions, lon, lat, scientificName, pasture, all_of(draw_col)) %>%
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
  forest_list[[paste0("forest_", draw, "_wide")]] <- dcast(as.data.table(forest_tmp), 
                                                           ecoregions + lon + lat ~ scientificName, 
                                                           value.var = draw_col, fill = 0)
  
  pasture_list[[paste0("pasture_", draw, "_wide")]] <- dcast(as.data.table(pasture_tmp), 
                                                             ecoregions + lon + lat ~ scientificName, 
                                                             value.var = draw_col, fill = 0)
}

# Construcción de grids
dggs <- lapply(resolutions, function(x) dgconstruct(res = x))
names(dggs) <- paste0("dggs_", resolutions)

# Convertir coordenadas a celdas
cells <- lapply(dggs, function(dg) dgGEO_to_SEQNUM(dg, points_coords$X, points_coords$Y)$seqnum)
uc <- lapply(cells, unique)

# Función para calcular perdidas locales y regionales de diversidad a través de iterciones posteriores

resultados <- list()
resultados_rasters <- list()
# Análisis por draw y resolución
for (draw in draws) {
  forest_data <- forest_list[[paste0("forest_", draw, "_wide")]]
  pasture_data <- pasture_list[[paste0("pasture_", draw, "_wide")]]
  
  # Inicializar estructura de resultados
  resultados[[paste0("draw_", draw)]] <- vector("list", length(resolutions))
  names(resultados[[paste0("draw_", draw)]]) <- paste0("res_", resolutions)
  
  cell_ratios <- get_avg_cell_ratios(forest_data, pasture_data, cutoff_type="absolute", cutoff=1)
  
  for (r in seq_along(resolutions)) {
    res <- resolutions[r]
    uc_res <- uc[[r]]
    cells_res <- cells[[r]]
    
    # Inicializar listas vacías para cada variable de interés
    cell_logratios <- numeric(length(uc_res))
    cell_richness <- numeric(length(uc_res))
    cell_pointwise_logratios <- numeric(length(uc_res))
    cell_pointwise_richness <- numeric(length(uc_res))
    cell_numbers <- numeric(length(uc_res))
    cell_beta <- numeric(length(uc_res))
    
    cell_logratios_expanded <- rep(NA, nrow(points_coords))
    cell_pointwise_logratios_expanded <- rep(NA, nrow(points_coords))
    cell_beta_expanded <- rep(NA, nrow(points_coords))
    
    
    for (i in seq_along(uc_res)) {
      cat("draw", draw, "res", res, "cell ", i, "\n")
      
      cell_positions <- which(cells_res == uc_res[i])
      grl <- get_regional_ratios(forest_data, pasture_data, cutoff_type="absolute", cutoff=1, cell_positions = cell_positions)
      
      cell_logratios[i] <- grl$med_logratio
      cell_richness[i] <- grl$n
      cell_logratios_expanded[which(cells_res == uc_res[i])] <- cell_logratios[i]
      cell_pointwise_logratios[i] <- mean(cell_ratios$med_logratio[which(cells_res == uc_res[i])], na.rm=TRUE)
      cell_pointwise_richness[i] <- mean(cell_ratios$n[cells_res == uc_res[i] & !is.na(cell_ratios$avg_logratio)])
      cell_beta[i] <- cell_richness[i] / cell_pointwise_richness[i]
      cell_beta_expanded[which(cells_res == uc_res[i])] <- cell_beta[i]
      cell_pointwise_logratios_expanded[which(cells_res == uc_res[i])] <- cell_pointwise_logratios[i]
      cell_numbers[i] <- sum(cells_res == uc_res[i])
    }
    
    # Crear rásteres
    raster_regional <- rasterFromXYZ(cbind(forest_data[, c("lon", "lat")], cell_logratios_expanded))
    raster_pointwise <- rasterFromXYZ(cbind(forest_data[, c("lon", "lat")], cell_pointwise_logratios_expanded))
    raster_beta <- rasterFromXYZ(cbind(forest_data[, c("lon", "lat")], cell_beta_expanded))
    raster_difference <- raster_regional - raster_pointwise
    
    # Guardar rásteres en la lista
    resultados_rasters[[paste0("draw_", draw)]][[paste0("res_", res)]] <- list(
      raster_regional = raster_regional,
      raster_pointwise = raster_pointwise,
      raster_beta = raster_beta,
      raster_difference = raster_difference
    )
  }
}


raster_to_df <- function(raster, resolution, draw, type) {
  df <- as.data.frame(raster, xy = TRUE)
  
  # Detectar la columna de logratio
  logratio_col <- grep("cell_.*_expanded|layer", colnames(df), value = TRUE)
  
  df %>%
    rename(logratio = all_of(logratio_col)) %>%
    mutate(res = paste0(resolution),
           draw = draw,
           type = type)
}

# Extraer datos de todos los draws y resoluciones
df_list <- list()

for (draw in names(resultados_rasters)) {
  for (res in names(resultados_rasters[[draw]])) {
    df_list[[paste0(draw, "_", res, "_regional")]] <-
      raster_to_df(resultados_rasters[[draw]][[res]]$raster_regional, resolution = res, draw = draw, type = "Regional")
    
    df_list[[paste0(draw, "_", res, "_pointwise")]] <-
      raster_to_df(resultados_rasters[[draw]][[res]]$raster_pointwise, resolution = res, draw = draw, type = "Pointwise")
    
    df_list[[paste0(draw, "_", res, "_difference")]] <-
      raster_to_df(resultados_rasters[[draw]][[res]]$raster_difference, resolution = res, draw = draw, type = "Difference")
    
    df_list[[paste0(draw, "_", res, "_beta")]] <-
      raster_to_df(resultados_rasters[[draw]][[res]]$raster_beta, resolution = res, draw = draw, type = "Beta")
  }
}

df_list <- lapply(df_list, function(df) {
  df$res <- as.character(df$res)  # Convertir a character
  return(df)
})

# Unir todo en un solo dataframe
df_all <- bind_rows(df_list) %>%
  filter(!is.na(logratio))

df_subset1 <- df_all %>% filter(type %in% c("Regional", "Pointwise"))
df_subset2 <- df_all %>% filter(type %in% c("Difference", "Beta"))

ggplot(df_subset1, aes(x = x, y = y, fill = logratio)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(type ~ res) + 
  theme_bw() +
  labs(title = "Comparison of regional and pointwise Log ratios", 
       fill = "Log ratio")


ggplot(df_subset2, aes(x = x, y = y, fill = logratio)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_grid(type ~ res, scales = "free") + 
  theme_bw() +
  labs(title = NULL, 
       fill = "ratio")

library(tidyr)

# Transformar el formato a "wide" para tener columnas separadas
df_subset2_wide <- df_subset2 %>%
  pivot_wider(names_from = type, values_from = logratio, 
              names_prefix = "") %>%
  rename(excess_loss = Difference, beta = Beta)

# Graficar la relación entre Beta y Excess Loss
ggplot(df_subset2_wide, aes(x = beta, y = excess_loss, group = draw)) +
  geom_hex(bins = 30) +  # Mapa de densidad hexagonal
  scale_fill_viridis_c(trans = "log") +  # Escala perceptual
  geom_smooth(method = "lm", color = "black", fill = "gray70", se = TRUE,  level = 0.95) +  # Línea de regresión con IC
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~res, nrow = 2) +  # Facetas por resolución
  theme_bw(base_size = 14) +
  labs(x = "Beta", y = "Excess regional loss", fill = "Density") +
  theme(legend.position = "none")

################################################################################
################################################################################

y# Plot the difference between the regional and the mean pointwise losses by cell
plot(raster_difference_10)
plot(raster_difference_7)plot(raster_difference_4)

# Check how the mean (across cells) of the regional losses decreases with decreasing cell size
mean(colombia_ratios_df$med_logratio)
mean(cell_logratios_4, na.rm=T)
mean(cell_logratios_6, na.rm=T)
mean(cell_logratios_7, na.rm=T)
mean(cell_logratios_8, na.rm=T)
mean(cell_logratios_9, na.rm=T)
mean(cell_logratios_10, na.rm=T)

mean(cell_ratios$med_logratio, na.rm = T)
mean(cell_pointwise_logratios_4)
mean(cell_pointwise_logratios_6)
mean(cell_pointwise_logratios_7)
mean(cell_pointwise_logratios_8)
mean(cell_pointwise_logratios_9)
mean(cell_pointwise_logratios_10)

library(tidyr)

df_corr <- df_all %>% 
  filter(type %in% c("Beta", "Difference")) %>% 
  pivot_wider(names_from = type, values_from = logratio)


ggplot(df_corr, aes(x = Beta, y = Difference)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Línea de regresión
  facet_wrap(~ res, scales = "free") +
  theme_bw() +
  labs(title = "Correlation between Beta and excess regonal loss",
       x = "Multiplicative beta-diversity",
       y = "Excess regonal loss")

# Comparación de pérdida regional vs. local para cada resolución
a <- (cell_logratios_10 - cell_pointwise_logratios_10)
b <- (cell_richness_10 / cell_pointwise_richness_10)

c <- (cell_logratios_9 - cell_pointwise_logratios_9)
d <- (cell_richness_9 / cell_pointwise_richness_9)

e <- (cell_logratios_8 - cell_pointwise_logratios_8)
f <- (cell_richness_8 / cell_pointwise_richness_8)

g <- (cell_logratios_7 - cell_pointwise_logratios_7)
h <- (cell_richness_7 / cell_pointwise_richness_7)

i <- (cell_logratios_6 - cell_pointwise_logratios_6)
j <- (cell_richness_6 / cell_pointwise_richness_6)

k <- (cell_logratios_4 - cell_pointwise_logratios_4)
l <- (cell_richness_4 / cell_pointwise_richness_4)

# Graficar las relaciones
plot(a ~ b, xlim = c(1, 8), ylim = c(-0.5, 1.5), main="Pérdida Regional vs. Relación de Riqueza", 
     xlab="Cociente Riqueza Regional/Local", ylab="Diferencia en Log-Ratio de Pérdida")
points(c ~ d, col = "blue")  # Resolución 9
points(e ~ f, col = "purple") # Resolución 8
points(g ~ h, col = "red")    # Resolución 7
points(i ~ j, col = "green")  # Resolución 6
points(k ~ l, col = "orange") # Resolución 4
points(6.5, .59, col = "black", pch = 16, cex = 1.2) # Punto de referencia

# Ajuste de modelos lineales para cada resolución
summary(lm(a ~ b))
summary(lm(c ~ d))
summary(lm(e ~ f))
summary(lm(g ~ h))
summary(lm(i ~ j))
summary(lm(k ~ l))



