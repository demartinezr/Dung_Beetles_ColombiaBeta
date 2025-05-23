setwd("C:/Users/PC/Dropbox/CO_DBdata")
library(sf)
library(dplyr)
library(tidyr)
library(tidyverse)
library(data.table)
library(purrr)
library(ggplot2)
library(ggridges)
library(patchwork)

ecoregions_predictions <- readRDS("./Analysis/mean_abundance/ecoregions_predictions_100.rds")

############# Analysis of species specific sensitivities #####################

################# local losers and winners

ecoregions_predictions <- map2(
  ecoregions_predictions, 
  names(ecoregions_predictions), 
  ~ mutate(.x, ecoregions = .y)
)

ecoregions_df <- bind_rows(ecoregions_predictions)
ecoregions_forest <- ecoregions_df %>% filter(pasture == 0) %>% select(lon, lat, ecoregions, scientificName, starts_with("abun__draw"), geometry)
ecoregions_pasture <- ecoregions_df %>% filter(pasture == 1) %>% select(lon, lat, ecoregions, scientificName, starts_with("abun__draw"), geometry)

cols_draw <- grep("^abun__draw", colnames(ecoregions_forest), value = TRUE)

ratio_values <- ecoregions_forest[, cols_draw] / ecoregions_pasture[, cols_draw]

#saveRDS(ratio_values, "./decline/ratio_values.rds")

colnames(ratio_values) <- gsub("^abun__", "ratio__", colnames(ratio_values))

ratio_col_df <- ecoregions_forest %>%
  select(lat, lon, ecoregions, scientificName) %>%
  bind_cols(ratio_values)

draw_summary_latlon <- ratio_col_df %>%
  group_by(lat_lon) %>%
  group_map(~ {
    df <- .x
    latlon <- df$lat_lon[1]
    geom <- df$geometry[1]  
    eco <- df$ecoregions[1]  
    
    result <- map_dfr(ratio_cols, function(col) {
      vals <- df[[col]]
      vals <- vals[is.finite(vals)]  
      
      tibble(
        draw = col,
        losers = sum(vals > 1),
        winners = sum(vals < 1),
        total = length(vals),
        lat_lon = latlon,
        ecoregions = eco
      )
    })
    # Añadir geometría al resumen usando el lat_lon
    result$geometry <- geom
    st_sf(result)
  }) %>%
  bind_rows()

# saveRDS(draw_summary_latlon, "./decline/draw_summary_latlon.rds")

draw_summary_latlon <- readRDS("./decline/draw_summary_latlon.rds")
draw_summary_latlon <- draw_summary_latlon %>% 
  select(-ecoregions) %>%
  st_drop_geometry(draw_summary_latlon)

draw_summary_latlon$ecoregions <- paste0("Local")

##############3#### ecoregions losers and winners

# Get the region names from the list names
eco_names <- names(ecoregions_predictions)
# calculate the proportional change in abundance for species and region by draw
abundance_change <- function(sf_df, eco_names) {
  sf_df %>%
    st_drop_geometry() %>% 
    # Divide into forest (pasture = 0) and pasture (pasture = 1)
    group_by(scientificName, pasture) %>%
    summarise(across(starts_with("abun__draw_"), \(x) sum(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_wider(names_from = pasture, values_from = starts_with("abun__draw_"), 
                names_glue = "{.value}_pasture{pasture}") %>%
    mutate(across(ends_with("_pasture1"), 
                  ~ get(sub("_pasture1", "_pasture0", cur_column())) / .x, 
                  .names = "ratio_{.col}")) %>%
    dplyr::select(scientificName, starts_with("ratio_")) %>%
    mutate(ecoregions = eco_names)
}
# Apply the function for each species
ratio_draw <- map2(ecoregions_predictions, eco_names, abundance_change)
# list to dataframe
ratio_draw_df <- bind_rows(ratio_draw) 

# loser and winners
ratio_cols <- grep("^ratio_abun__.*_pasture1$", names(ratio_draw_df), value = TRUE)

# Agrupar por ecoregions y contar por draw
draw_summary_ecoregions <- ratio_draw_df %>%
  group_by(ecoregions) %>%
  group_map(~ {
    df <- .x
    eco <- .y$ecoregions
    map_dfr(ratio_cols, function(col) {
      vals <- df[[col]]
      vals <- vals[is.finite(vals)]  # Eliminar NaN, Inf, 0, negativos
      
      tibble(
        draw = col,
        losers = sum(vals > 1),
        winners = sum(vals < 1),
        total = length(vals),
        ecoregions = eco
      )
    })
  }) %>% bind_rows()

#saveRDS(draw_summary_ecoregions, "./decline/draw_summary_ecoregions.rds")

#draw_plot_df <- draw_summary_ecoregions %>%
#  mutate(
#    pct_winners = winners / total * 100,
#    pct_losers = losers / total * 100
#  ) %>%
#  select(ecoregions, pct_winners, pct_losers) %>%
#  gather(key = "group", value = "percentage", pct_winners, pct_losers)
draw_summary_ecoregions <- readRDS("./decline/draw_summary_ecoregions.rds")

draw_summary_ecoregions$ecoregions <- gsub(" forests", "", draw_summary_ecoregions$ecoregions)
draw_summary_ecoregions$ecoregions <- recode(draw_summary_ecoregions$ecoregions,
                                  "Apure-Villavicencio dry"         = "Villavicencio dry",
                                  "Cauca Valley montane"            = "Cauca montane",
                                  "Cordillera Oriental montane"     = "EC montane",
                                  "Eastern Cordillera Real montane" = "CC montane",
                                  "Magdalena Valley dry"            = "Magdalena dry",
                                  "Magdalena Valley montane"        = "Magdalena Valley",
                                  "Northern Andean páramo"          = "Andean páramo",
                                  "Northwest Andean montane"        = "WC montane")

draw_summary_ecoregions$draw <- gsub("ratio_abun__draw_(\\d+)_pasture1", "ratio__draw_\\1", draw_summary_ecoregions$draw)

#draw_plot_df <-  draw_summary_ecoregions %>%
#  pivot_longer(cols = c(losers, winners),
#               names_to = "category",
#               values_to = "n_species")

#ecoregion_order <- draw_summary_ecoregions %>%
#  group_by(ecoregions) %>%
#  summarise(mean_losers = mean(losers)) %>%
#  arrange(mean_losers) %>%
#  pull(ecoregions)

#draw_plot_df <- draw_plot_df %>%
#  mutate(ecoregions = factor(ecoregions, levels = ecoregion_order))


    # Graficar
#ggplot(draw_plot_df, aes(x = ecoregions, y = n_species, fill = category)) +
#  geom_violin(trim = FALSE, alpha = 0.7) +
#  geom_jitter(aes(color = category),
#              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
#              size = 1.5, alpha = 0.5) +
#  labs(x = "Ecorregión", y = "Número de especies", fill = "Categoría", color = "Categoría") +
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ecoregion_avg <- draw_summary_ecoregions %>%
  group_by(ecoregions) %>%
  summarise(mean_losers = mean(losers),
            mean_winners = mean(winners)) %>%
  mutate(total_mean = mean_losers + mean_winners,
         perc_losers = mean_losers / total_mean,
         perc_winners = mean_winners / total_mean) %>%
  select(ecoregions, perc_losers, perc_winners)

# Pasar a formato largo
ordered_levels <- ecoregion_avg$ecoregions[order(-ecoregion_avg$perc_losers)]
ecoregion_avg$ecoregions <- factor(ecoregion_avg$ecoregions, levels = ordered_levels)

ecoregion_avg_long <- ecoregion_avg %>%
  pivot_longer(cols = starts_with("perc_"),
               names_to = "category",
               values_to = "proportion") %>%
  mutate(category = recode(category,
                           perc_losers = "Losers",
                           perc_winners = "Winners"))


# Graficar barras apiladas
fig_2c <- ggplot(ecoregion_avg_long, aes(x = proportion, y = ecoregions, fill = category)) +
              geom_bar(stat = "identity") +
              scale_fill_viridis_d(direction = -1) +
              labs(title ="c", x = "Proportion of species (Regional)", y = NULL, fill = "Group") +
              theme_bw() +
  theme(plot.title = element_text(size = 13),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        axis.title.x = element_text(size = 13),
        legend.position = "none")

### local ecoregion losers and winners
draw_summary_latlon_eco <- readRDS("./decline/draw_summary_latlon.rds")

ecoregion_pointwise <- draw_summary_latlon_eco %>%
  group_by(ecoregions) %>%
  summarise(mean_losers = mean(losers),
            mean_winners = mean(winners)) %>%
  mutate(total_mean = mean_losers + mean_winners,
         perc_losers = mean_losers / total_mean,
         perc_winners = mean_winners / total_mean) %>%
  select(ecoregions, perc_losers, perc_winners)
#saveRDS(ecoregion_pointwise, "./decline/ecoregion_pointwise.rds")

ecoregion_pointwise <- readRDS("./decline/ecoregion_pointwise.rds")
ecoregion_pointwise$ecoregions <- gsub(" forests", "", ecoregion_pointwise$ecoregions)
ecoregion_pointwise$ecoregions <- recode(ecoregion_pointwise$ecoregions,
                                             "Apure-Villavicencio dry"         = "Villavicencio dry",
                                             "Cauca Valley montane"            = "Cauca montane",
                                             "Cordillera Oriental montane"     = "EC montane",
                                             "Eastern Cordillera Real montane" = "CC montane",
                                             "Magdalena Valley dry"            = "Magdalena dry",
                                             "Magdalena Valley montane"        = "Magdalena Valley",
                                             "Northern Andean páramo"          = "Andean páramo",
                                             "Northwest Andean montane"        = "WC montane")
# Pasar a formato largo
ordered_levels <- ecoregion_pointwise$ecoregions[order(-ecoregion_pointwise$perc_losers)]
ecoregion_pointwise$ecoregions <- factor(ecoregion_pointwise$ecoregions, levels = ordered_levels)

ecoregion_pointwise_long <- ecoregion_pointwise %>%
  pivot_longer(cols = starts_with("perc_"),
               names_to = "category",
               values_to = "proportion") %>%
  mutate(category = recode(category,
                           perc_losers = "Losers",
                           perc_winners = "Winners"))

# Graficar barras apiladas
fig_2b <- ggplot(ecoregion_pointwise_long, aes(x = proportion, y = ecoregions, fill = category)) +
              geom_bar(stat = "identity") +
              scale_fill_viridis_d(direction = -1) +
              labs(title="b", x = expression("Proportion of species (2 km"^2*")"), y = NULL, fill = "Group") +
              theme_bw() +
              theme(plot.title = element_text(size = 13),
                    axis.text.x = element_text(size = 12, color = "black"),
                    axis.text.y = element_text(size = 12, color = "black"), 
                    axis.title.x = element_text(size = 13),
                    legend.position = "none")

##### Pan Colombia

species_means <- ratio_draw_df %>%
  group_by(scientificName) %>%
  summarise(across(all_of(ratio_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Calcular resumen por especie (winners, losers, total y promedio del ratio)
draw_summary_sp <- map_dfr(ratio_cols, function(col) {
  vals <- species_means[[col]]
  vals <- vals[is.finite(vals)]  # NaN, Inf
  
  tibble(
    draw = col,
    losers = sum(vals > 1),
    winners = sum(vals < 1),
    total = length(vals),
    ecoregions="Pan-Colombia"
#    prop_losers = sum(vals > 1) / length(vals),
#    prop_winners = sum(vals < 1) / length(vals)
  )
})

#saveRDS(draw_summary_sp, "./decline/draw_summary_sp.rds")

draw_summary_sp <- readRDS("./decline/draw_summary_sp.rds")
draw_summary_sp$draw <- gsub("ratio_abun__draw_(\\d+)_pasture1", "ratio__draw_\\1", draw_summary_sp$draw)

draw_summary_all <- bind_rows(draw_summary_sp, draw_summary_latlon)

col_avg <- draw_summary_all %>%
  group_by(ecoregions) %>%
  summarise(mean_losers = mean(losers),
            mean_winners = mean(winners)) %>%
  mutate(total_mean = mean_losers + mean_winners,
         perc_losers = mean_losers / total_mean,
         perc_winners = mean_winners / total_mean) %>%
  select(ecoregions, perc_losers, perc_winners)

# Pasar a formato largo
ordered_levels_col <- col_avg$ecoregions[order(-col_avg$perc_losers)]
col_avg$ecoregions <- factor(col_avg$ecoregions, levels = ordered_levels_col)

col_avg_long <- col_avg %>%
  pivot_longer(cols = starts_with("perc_"),
               names_to = "category",
               values_to = "percentage") %>%
  mutate(category = recode(category,
                           perc_losers = "Losers",
                           perc_winners = "Winners"))

# Graficar barras apiladas
fig_2a <- ggplot(col_avg_long, aes(x = percentage, y = ecoregions, fill = category)) +
              geom_bar(stat = "identity", alpha=0.9) +
              scale_fill_viridis_d(direction = -1) +
              labs(title="a", x = "Proportion of species", y = NULL, fill = NULL) +
              theme_bw() +
              theme(plot.title = element_text(size = 13),
                    axis.text.x = element_text(size = 12, color = "black"),
                    axis.text.y = element_text(size = 12, color = "black"), 
                    axis.title.x = element_text(size = 13),
                    legend.text = element_text(size = 13),
                    legend.key.size = unit(1.5, "lines"))


fila_1 <- plot_spacer() + fig_2a + plot_spacer()
fila_1 <- fila_1 + plot_layout(ncol = 3, widths = c(0.5, 2, 0.5))
fila_2 <- fig_2b + fig_2c
fig_2 <- fila_1 / fila_2 + plot_layout(heights = c(0.2, 0.6))

ggsave("./fig_2.jpeg", plot = fig_2, width = 10, height = 5.5, units = "in",        
       dpi = 300, device = "jpeg")

fig_2_2 <- fig_2a / fig_2b / fig_2c + plot_layout(heights = c(0.8, 3, 3))
ggsave("./fig_2_2.jpeg", plot = fig_2_2, width = 8.5, height = 11, units = "in",        
       dpi = 300, device = "jpeg")
##############################################################################
ecoregions_all <- as.data.frame(readRDS("./ecoregions_100.rds"))

# function to compute biodiversity loss  
source("F:/repositorio/Dung_Beetles_ColombiaBeta/9_biodiversity_loss/db_compute_loss.R")
draws <-  30 * c(1:100) # 100 draws
unique_ecoregions <- unique(ecoregions_all$ecoregions)

colombia_decline <- list()
local_decline <- list()
ecoregion_decline <- list()

for (i in seq_along(draws)) {
  cat("draw:", i, "of", length(draws), "\n")  
  draw_col <- paste0("abun__draw_", draws[i])
  
  # Preparar los datos
  data_draw <- ecoregions_all %>%
    select(ecoregions, lon, lat, scientificName, pasture, all_of(draw_col)) %>%
    mutate(cell_id = row_number()) %>%
    select(cell_id, everything())
  row.names(data_draw) <- NULL
  
  # Dividir en bosque y pasto
  forest_draw <- data_draw %>% filter(pasture == 0) %>% select(-pasture)
  pasture_draw <- data_draw %>% filter(pasture == 1) %>% select(-pasture)
  
  # Convertir a formato ancho
  forest_draw_dt <- as.data.table(forest_draw)
  forest_draw_wide <- dcast(forest_draw_dt, ecoregions + lon + lat ~ scientificName, 
                            value.var = draw_col, fill = 0)
  
  pasture_draw_dt <- as.data.table(pasture_draw)
  pasture_draw_wide <- dcast(pasture_draw_dt, ecoregions + lon + lat ~ scientificName, 
                             value.var = draw_col, fill = 0)
  
  # 1. Declive nacional (Colombia completo)
  colombia_result <- get_sample_percent_decline(forest_draw_wide, pasture_draw_wide, cutoff = 1)
  colombia_decline[[draw_col]] <- colombia_result
  
  # 2. Declive local (celda por celda)
  cell_positions_all <- 1:nrow(forest_draw_wide)
  local_result <- get_sample_percent_decline(forest_draw_wide, pasture_draw_wide, cutoff = 1, cell_positions = cell_positions_all)
  local_decline[[draw_col]] <- local_result$pointwise
  
  # 3. Declive por ecorregión
  ecoregion_results <- list()
  for (j in seq_along(unique_ecoregions)) {
    cat("  Ecoregion:", j, "of", length(unique_ecoregions), "\n")
    
    eco <- unique_ecoregions[j]
    cell_positions <- which(forest_draw_wide$ecoregions == eco)
    
    eco_result <- get_sample_percent_decline(forest_draw_wide, pasture_draw_wide, 
                                             cutoff = 1, 
                                             cell_positions = cell_positions)
    eco_result$ecoregion <- eco
    ecoregion_results[[eco]] <- eco_result
  }
  ecoregion_decline[[draw_col]] <- ecoregion_results
  
  # Liberar memoria
  rm(forest_draw, forest_draw_dt, forest_draw_wide, 
     pasture_draw, pasture_draw_dt, pasture_draw_wide, data_draw)
  gc()
}

saveRDS(colombia_decline, "./decline/colombia_decline.rds")
saveRDS(local_decline, "./decline/local_decline.rds")
saveRDS(ecoregion_decline, "./decline/ecoregion_decline.rds")

###############################################################################
colombia_decline <- readRDS("./decline/colombia_decline.rds")
local_decline <- readRDS("./decline/local_decline.rds")
ecoregion_decline <- readRDS("./decline/ecoregion_decline.rds")

# regional and local mean decline
mean_decline <- do.call(rbind, lapply(seq_along(colombia_decline), function(i) {
  draw <- colombia_decline[[i]]
  data.frame(
    draw = names(colombia_decline)[i],
    ecoregion = "Colombia",
    Local = mean(draw$pointwise, na.rm = TRUE),
    Regional = draw$total,
    nsp = draw$nsp
  )
}))


# Convertimos a formato largo para ggplot
mean_long <- tidyr::pivot_longer(mean_decline, cols = c(Local, Regional),
                                 names_to = "scale", values_to = "decline")

# Gráfico de violín
ggplot(mean_long, aes(x = scale, y = decline, fill = scale)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +
  xlab("") +
  ylab("Proportion of species with decline > 1") +
  theme_bw() +
  theme(legend.position = "none")

mean_decline$difference <- mean_decline$Regional / mean_decline$Local

ggplot(mean_decline, aes(x = "", y = difference)) +
  geom_violin(fill = "steelblue", alpha = 0.8, color = NA) +
  geom_jitter(width = 0.05, alpha = 0.5, size = 1) +
  ylab("Regional - Local decline") +
  xlab("") +
  theme_bw(base_size = 14)
#############################################################################

ecoregion_mean_decline <- rbindlist(
  lapply(names(ecoregion_decline), function(draw_name) {
    draw <- ecoregion_decline[[draw_name]]
    data.table(
      draw = draw_name,
      ecoregion = names(draw),
      Local = sapply(draw, function(x) mean(x$pointwise, na.rm = TRUE)),
      Regional = sapply(draw, function(x) x$total),
      nsp = sapply(draw, function(x) x$nsp)
    )
  })
)
ecoregion_mean_decline$ecoregion <- gsub(" forests", "", ecoregion_mean_decline$ecoregion)
ecoregion_mean_decline$ecoregion <- recode(ecoregion_mean_decline$ecoregion,
                         "Apure-Villavicencio dry"         = "Villavicencio dry",
                         "Cauca Valley montane"            = "Cauca montane",
                         "Cordillera Oriental montane"     = "EC montane",
                         "Eastern Cordillera Real montane" = "CC montane",
                         "Magdalena Valley dry"            = "Magdalena dry",
                         "Magdalena Valley montane"        = "Magdalena Valley",
                         "Northern Andean páramo"          = "Andean páramo",
                         "Northwest Andean montane"        = "WC montane")

ecoregion_long <- melt(ecoregion_mean_decline,
                       id.vars = c("draw", "ecoregion", "nsp"),
                       measure.vars = c("Local", "Regional"),
                       variable.name = "scale",
                       value.name = "decline")


ecoregion_long$scale <- factor(ecoregion_long$scale, levels = c("Local", "Regional"))

ecoregion_mean_decline$difference <- ecoregion_mean_decline$Regional / ecoregion_mean_decline$Local


# join pan colombia and ecoregions

decline_all <- rbind(ecoregion_long, mean_long)

regional_order <- decline_all %>%
  filter(scale == "Regional") %>%
  arrange(desc(decline)) %>%
  pull(ecoregion)

ecoregion_levels <- unique(c(regional_order, 
                             decline_all %>% filter(scale == "Local") %>% pull(ecoregion)))

decline_all <- decline_all %>%
  mutate(ecoregion = factor(ecoregion, levels = ecoregion_levels))


#fig_2b <- ggplot(decline_all, aes(x = ecoregion, y = decline, fill = scale)) +
#              geom_violin(trim = FALSE, alpha = 0.4, width = 0.9, position = position_dodge(width = 0.8), linewidth = 0.5) +
#              geom_boxplot(width = 0.2, color = "black", alpha = 0.4, outlier.shape = NA, position = position_dodge(width = 0.8), linewidth = 0.5) +
#              scale_fill_viridis_d() +           
#              theme_bw() + 
#              scale_y_continuous(limits=c(0.4, 1)) +
#              theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, color="black"),
#                    axis.text.y = element_text(size = 10, color="black"), 
#                    axis.title.x = element_text(size = 11),
#                    axis.title.y = element_text(size = 11),
#                    legend.position = "bottom",
#                    legend.text = element_text(size = 11),
#                    legend.title = element_text(size = 12)) +
#              labs(x = NULL,
#                   y = "Proportion of declining species") +
#              geom_vline(xintercept = seq(1, length(ecoregion_levels) - 1) + 0.5, 
#                         color = "gray40", linetype = "solid", linewidth = 0.4) + coord_flip()


fig_3a <- ggplot(decline_all, aes(x = decline, y = reorder(ecoregion, -decline), fill = scale)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5)  +
  scale_fill_viridis_d() +
  labs(title = "a", x = "Proportion of declining species", y = NULL) +
  theme_bw() +
  theme(legend.position = "right",
        plot.title = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(size = 11, color = "black"),
              axis.text.y = element_text(size = 11, color = "black"), 
              axis.title.x = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 13)) +
  scale_x_continuous(limits=c(0.5, 1)) +
  stat_density(geom = "line", position = "identity", aes(x = decline), kernel = "gaussian")


diff_all <- rbind(mean_decline, ecoregion_mean_decline)

#fig_2c <- ggplot(diff_all, aes(reorder(x = ecoregion, difference), y = difference, fill = ecoregion)) +
#  geom_violin(trim = FALSE, alpha = 0.8, width = 0.8, position = position_dodge(width = 0.8)) +
#  geom_boxplot(width = 0.2, color = "black", alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
#  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.7) +
#  scale_fill_viridis_d() +
#  theme_bw() + 
#  scale_y_continuous(limits=c(0.7, 2)) +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, color="black"),
#        axis.text.y = element_text(size = 10, color="black"),
#        axis.title.y = element_text(size = 11),
#        legend.position = "none") +
#  labs(x = NULL,
#       y = "Proportion of declining species")

fig_3b <- ggplot(diff_all, aes(x = difference, y = reorder(ecoregion, difference), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5)  +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
  scale_fill_viridis_d() +
  labs(title = "b", x = "Declining species proportion ratio\n(Regional/local)", y = NULL) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title.x = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 12, face = "bold")) +
  guides(fill = "none") + 
  scale_x_continuous(limits=c(0.75, 1.75)) +
  stat_density(geom = "line", position = "identity", aes(x = difference), kernel = "gaussian")

fig_3_2 <- fig_3a / fig_3b


ggsave("./fig_3_2.jpeg", plot = fig_3_2, width = 6.5, height = 8.5, units = "in",        
       dpi = 300, device = "jpeg")


decline_summary <- decline_all %>%
  group_by(ecoregion, scale) %>%
  summarise(
    mean_decline = mean(decline, na.rm = TRUE),
    lower_ci = quantile(decline, 0.025, na.rm = TRUE),
    upper_ci = quantile(decline, 0.975, na.rm = TRUE),
    .groups = "drop"
  )




