# ------------------------------------------------------------------------------
# This script implements the richness-based biodiversity loss metric used in the 
# manuscript to quantify the proportion of species declining following forest–
# pasture conversion.
#
# Using posterior predictions from the hierarchical abundance model, the script
# compares predicted species abundances in forest and pasture to calculate
# species-specific sensitivity ratios (forest / pasture). Species with ratios
# greater than 1 are classified as declining (“losers”), indicating lower
# predicted abundance in pasture.The proportion of declining species is summarized 
# across three spatial scales
#
# The workflow is divided into two main analyses:
#
#   1) Richness-based biodiversity loss across spatial scales (Fig. 2a)
#   2) Spatial scaling of richness-based biodiversity loss across ecoregions (Fig. 3)
# ------------------------------------------------------------------------------
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
library(RColorBrewer)

# Load abundance predictions dataset
ecoregions_predictions <- readRDS("./Analysis/mean_abundance/ecoregions_predictions_100_mod.rds")

# ------------------------------------------------------------------------------
# Compute species-level forest–pasture abundance ratios across posterior draws
#
# Predicted abundances for each species are available for forest and pasture
# scenarios across grid cells and posterior draws. Here we calculate the
# forest-to-pasture abundance ratio for each species and draw, which represents
# species-level sensitivity to land-use conversion.

# Add ecoregion identifiers to prediction data frames and combine them
ecoregions_predictions <- map2(
  ecoregions_predictions,names(ecoregions_predictions), 
  ~ mutate(.x, ecoregions = .y))

# Separate predicted abundances for forest and pasture scenarios
ecoregions_df <- bind_rows(ecoregions_predictions)

ecoregions_forest <- ecoregions_df %>% filter(pasture == 0) %>% 
  select(lon, lat, ecoregions, scientificName, starts_with("abun__draw"), geometry)
ecoregions_pasture <- ecoregions_df %>% filter(pasture == 1) %>% 
  select(lon, lat, ecoregions, scientificName, starts_with("abun__draw"), geometry)

# Identify posterior draw columns
cols_draw <- grep("^abun__draw", colnames(ecoregions_forest), value = TRUE)

ratio_values <- ecoregions_forest[, cols_draw] / ecoregions_pasture[, cols_draw]
#saveRDS(ratio_values, "./decline/ratio_values_mod.rds")

# Rename columns to reflect ratio values instead of predicted abundances
colnames(ratio_values) <- gsub("^abun__", "ratio__", colnames(ratio_values))
# Combine ratio values with spatial and species identifiers
ratio_col <- ecoregions_forest %>%
  select(lat, lon, ecoregions, scientificName) %>%
  bind_cols(ratio_values)

#saveRDS(ratio_col_df, "./decline/ratio_col_df_mod.rds")

# Load ratio data
ratio_col_df <- readRDS("./decline/ratio_col_df_mod.rds")

# Identify ratio columns corresponding to posterior draws
ratio_cols <- grep("^ratio__draw_", colnames(ratio_col_df), value = TRUE)

# Create unique identifiers for each grid cell
ratio_col_df$lat_lon <- paste0(ratio_col_df$lat, "_", ratio_col_df$lon)

# Convert to spatial object
ratio_col_df <- st_as_sf(ratio_col_df, coords = c("lon", "lat"), crs = 4326)

# ------------------------------------------------------------------------------
# Summarize species responses within grid cells across posterior draws
#
# For each grid cell and posterior draw, species are classified as:
#   - losers: species with forest–pasture ratio > 1
#   - winners: species with ratio < 1
#
# The number of losers and winners species is counted to generate the
# richness-based biodiversity loss metric (proportion of declining species).

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
    result$geometry <- geom
    st_sf(result)
  }) %>%
  bind_rows()

#saveRDS(draw_summary_latlon, "./decline/draw_summary_latlon_eco_mod.rds")

################################################################################
# 1) Compute richness-based biodiversity loss metrics across spatial scales
#
# Using the counts of losing and winning species for each posterior draw,
# we calculate the proportion of declining species (losers / total species)
# within spatial units. These values are then aggregated to obtain mean
# community responses at three spatial scales:
#
#   1. Local grid cells (2 × 2 km)
#   2. Ecoregions
#   3. Near-national (pan-Colombia)
#
# These summaries form the basis of the richness-based biodiversity loss
# metric presented in Figure 2a.

draw_summary_latlon <- readRDS("./decline/draw_summary_latlon_eco_mod.rds")

# ------------------------------------------------------------------------------
# Calculate the proportion of declining and increasing species within
# each grid cell and posterior draw, then average across grid cells
# to obtain local-scale estimates.

draw_summary_latlon_mod <- draw_summary_latlon %>%
  st_drop_geometry() %>% 
  mutate(prop_losers = losers / (losers + winners),
    prop_winners = winners / (losers + winners)) %>%
  group_by(draw) %>%
  summarise(
    prop_losers = mean(prop_losers, na.rm = TRUE),
    prop_winners = mean(prop_winners, na.rm = TRUE))

# ------------------------------------------------------------------------------
# Ecoregion-scale biodiversity loss
#
# Species responses are aggregated across grid cells within each ecoregion.
# The proportion of declining species is calculated for each posterior draw
# and averaged across ecoregions.

draw_summary_ecoregions <- readRDS("./decline/draw_summary_ecoregions_eco_mod.rds")
draw_summary_ecoregions_mod <- draw_summary_ecoregions %>%
  mutate(prop_losers = losers / (losers + winners),
         prop_winners = winners / (losers + winners))  %>%
  group_by(draw) %>%
  summarise(
    prop_losers = mean(prop_losers, na.rm = TRUE),
    prop_winners = mean(prop_winners, na.rm = TRUE))

draw_summary_ecoregions_mod$draw <- gsub("ratio_abun__draw_(\\d+)_pasture1", "ratio__draw_\\1", 
                                         draw_summary_ecoregions_mod$draw)
draw_summary_ecoregions$ecoregions <- gsub(" forests", "", draw_summary_ecoregions$ecoregions)

# Standardize ecoregion names for consistency in plots and summaries
draw_summary_ecoregions$ecoregions <- 
  recode(draw_summary_ecoregions$ecoregions,
         "Apure-Villavicencio dry"         = "Villavicencio dry",
         "Cauca Valley montane"            = "Cauca montane",
         "Caqueta moist"                   ="Caquetá moist",   
         "Cordillera Oriental montane"     = "EC montane",
         "Eastern Cordillera Real montane" = "CC montane",
         "Magdalena Valley dry"            = "Magdalena dry",
         "Magdalena Valley montane"        = "Magdalena montane",
         "Northern Andean páramo"          = "Andean páramo",
         "Northwest Andean montane"        = "WC montane")

# ------------------------------------------------------------------------------
# Near-national biodiversity loss
#
# Predicted abundances are pooled across all grid cells to estimate the
# overall proportion of declining species across the modeled extent
# (near-national scale).

draw_summary_pc <- readRDS("./decline/draw_summary_sp.rds")
draw_summary_pc$draw <- gsub("ratio_abun__draw_(\\d+)_pasture1", "ratio__draw_\\1", draw_summary_pc$draw)
draw_summary_pc_mod <- draw_summary_pc %>%
  mutate(prop_losers = losers / (losers + winners),
         prop_winners = winners / (losers + winners))
draw_summary_pc_mod$ecoregions <- gsub(" forests", "", draw_summary_pc_mod$ecoregions)
draw_summary_pc_mod$ecoregions <- 
  recode(draw_summary_pc_mod$ecoregions,
         "Apure-Villavicencio dry"         = "Villavicencio dry",
         "Cauca Valley montane"            = "Cauca montane",
         "Caqueta moist"                   ="Caquetá moist",   
         "Cordillera Oriental montane"     = "EC montane",
         "Eastern Cordillera Real montane" = "CC montane",
         "Magdalena Valley dry"            = "Magdalena dry",
         "Magdalena Valley montane"        = "Magdalena montane",
         "Northern Andean páramo"          = "Andean páramo",
         "Northwest Andean montane"        = "WC montane")
# ------------------------------------------------------------------------------
# Calculate mean biodiversity loss estimates reported in the Results section

# local scale 2km metrics
mean(draw_summary_latlon_mod$prop_losers, na.rm = T)
# ecoregions metrics
mean(draw_summary_ecoregions_mod$prop_losers, na.rm = T)
# Pan-Colombia metrics  
mean(draw_summary_pc_mod$prop_losers, na.rm = T)

# ------------------------------------------------------------------------------
# Combine biodiversity loss estimates across spatial scales
#
# Posterior draws from each spatial scale are combined into a single
# dataset used to visualize the distribution of richness-based biodiversity
# loss across scales (Figure 2a).

richness_scales <- bind_rows(
  draw_summary_latlon_mod %>%
    group_by(draw) %>%
    summarise(across(c(prop_losers, prop_winners), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop") %>%
    mutate(scale = "Local"),
  
draw_summary_ecoregions_mod %>%
    mutate(scale = "Ecoregion"),
  
draw_summary_pc_mod %>%
    mutate(scale = "Near national"))

# ------------------------------------------------------------------------------
# Figure 2a: Richness-based biodiversity loss across spatial scales
#
# Violin plots show posterior distributions of the proportion of declining
# species at local, ecoregional, and near-national scales. Grey points
# represent individual posterior draws.

richness_scales$scale <- factor(richness_scales$scale,
                                levels = c("Local", "Ecoregion", "Near national"))


colors_scales <- c("Local"= "#762a83", "Ecoregion" = "#1b7837", "Near national" = "yellow")

fig_2a <- ggplot(richness_scales, aes(x = reorder(scale, prop_losers), y = prop_losers, fill = scale)) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.4, color="gray80") +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_fill_manual(values = colors_scales) +
  scale_y_continuous(limits=c(0.2, 1), breaks = seq(0.2, 1, by = 0.2)) +
  labs(title = "a", x = NULL, y = "Proportion of declining species", fill = "Scale") +
  theme_bw() +
  theme(legend.position = "none",
        legend.position.inside  = c(0.99, 0.3),         
        legend.justification = c("right", "bottom"),  
        legend.direction = "horizontal",
        legend.box.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
        legend.box.margin = margin(2, 2, 2, 2),
        plot.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 13, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_blank())

fig_2 <- fig_2a / fig_2b + plot_layout(heights = c(1, 1))

ggsave("./fig_2.pdf", plot = fig_2, width = 6, height = 7, units = "in",  dpi = 300)

################################################################################
# 2) Spatial scaling of richness-based biodiversity loss across ecoregions
#
# This analysis quantifies how the proportion of species declining under
# forest–pasture conversion changes with spatial scale across ecoregions.
#
# Using posterior predictions from the hierarchical abundance model, species-
# level forest–pasture abundance ratios are used to classify species as
# declining when predicted abundance is higher in forest than in pasture.
#
# For each posterior draw, the proportion of declining species is calculated
# at three spatial grains:
#
#   - Local grid cells (2 × 2 km)
#   - Ecoregional assemblages
#   - The near-national extent (all grid cells combined)
#
# These estimates allow evaluation of whether biodiversity loss inferred
# from local communities underestimates regional biodiversity responses
# to land-use conversion.
#
# The resulting posterior distributions are used to generate the density
# plots comparing local and ecoregional biodiversity loss across ecoregions
# shown in Figure 3.
#
# ------------------------------------------------------------------------------
# Load predicted species abundances for forest and pasture scenarios across
# grid cells and posterior draws
ecoregions_all <- as.data.frame(readRDS("./ecoregions_100.rds"))

# Load functions used to compute richness-based biodiversity loss metrics
source("F:/repositorio/Dung_Beetles_ColombiaBeta/9_biodiversity_loss/db_compute_loss.R")

# Select posterior draws from the abundance model used in the analysis
draws <-  30 * c(1:100)

# Identify the set of ecoregions represented in the predictions
unique_ecoregions <- unique(ecoregions_all$ecoregions)

# Initialize containers to store biodiversity loss estimates
# at near-national, local, and ecoregional spatial scales
colombia_decline <- list()
local_decline <- list()
ecoregion_decline <- list()

# Iterate across posterior draws to propagate model uncertainty
for (i in seq_along(draws)) {
  cat("draw:", i, "of", length(draws), "\n")  
  draw_col <- paste0("abun__draw_", draws[i])
  
  # Extract predicted abundances corresponding to the current posterior draw
  # and retain spatial and species identifiers
  data_draw <- ecoregions_all %>%
    select(ecoregions, lon, lat, scientificName, pasture, all_of(draw_col)) %>%
    mutate(cell_id = row_number()) %>%
    select(cell_id, everything())
  row.names(data_draw) <- NULL
  
  # Separate predicted abundances for forest and pasture scenarios
  forest_draw <- data_draw %>% filter(pasture == 0) %>% select(-pasture)
  pasture_draw <- data_draw %>% filter(pasture == 1) %>% select(-pasture)
  
  # Convert long-format predictions into species-by-site community matrices
  # required for biodiversity loss calculations
  forest_draw_dt <- as.data.table(forest_draw)
  forest_draw_wide <- dcast(forest_draw_dt, ecoregions + lon + lat ~ scientificName, 
                            value.var = draw_col, fill = 0)
  
  pasture_draw_dt <- as.data.table(pasture_draw)
  pasture_draw_wide <- dcast(pasture_draw_dt, ecoregions + lon + lat ~ scientificName, 
                             value.var = draw_col, fill = 0)
  
  # ---------------------------------------------------------------------------
  # Near-national biodiversity loss
  #
  # Calculate the proportion of species declining (forest > pasture)
  # across the entire modeled extent.

  colombia_result <- get_sample_percent_decline(forest_draw_wide, pasture_draw_wide, cutoff = 1)
  colombia_decline[[draw_col]] <- colombia_result
  
  # ---------------------------------------------------------------------------
  # Local biodiversity loss (2 × 2 km grid cells)
  #
  # Calculate the proportion of declining species within each grid cell.
  cell_positions_all <- 1:nrow(forest_draw_wide)
  local_result <- get_sample_percent_decline(forest_draw_wide, pasture_draw_wide,
                                             cutoff = 1, cell_positions = cell_positions_all)
  local_decline[[draw_col]] <- local_result$pointwise
  
  # ---------------------------------------------------------------------------
  # Ecoregional biodiversity loss
  #
  # For each ecoregion, aggregate the grid cells belonging to that region
  # and compute the proportion of species declining under forest–pasture
  # conversion.
  ecoregion_results <- list()
  for (j in seq_along(unique_ecoregions)) {
    cat("  Ecoregion:", j, "of", length(unique_ecoregions), "\n")
    eco <- unique_ecoregions[j]
    
    # Identify grid cells belonging to the current ecoregion
    cell_positions <- which(forest_draw_wide$ecoregions == eco)
    
    eco_result <- get_sample_percent_decline(forest_draw_wide, pasture_draw_wide, 
                                             cutoff = 1, 
                                             cell_positions = cell_positions)
    eco_result$ecoregion <- eco
    ecoregion_results[[eco]] <- eco_result}
  ecoregion_decline[[draw_col]] <- ecoregion_results
  
  # Remove temporary objects to reduce memory usage during iteration
  rm(forest_draw, forest_draw_dt, forest_draw_wide, 
     pasture_draw, pasture_draw_dt, pasture_draw_wide, data_draw)
  gc()}

# Save biodiversity loss estimates for subsequent analyses and visualization
saveRDS(colombia_decline, "./decline/colombia_decline_eco.rds")
saveRDS(local_decline, "./decline/local_decline_eco.rds")
saveRDS(ecoregion_decline, "./decline/ecoregion_decline_eco.rds")
#
# ------------------------------------------------------------------------------
# Load previously computed decline estimates for each spatial scale
colombia_decline <- readRDS("./decline/colombia_decline_eco.rds")
local_decline <- readRDS("./decline/local_decline_eco.rds")
ecoregion_decline <- readRDS("./decline/ecoregion_decline_eco.rds")
#
# ------------------------------------------------------------------
# Compute mean biodiversity loss across posterior draws
#
# For each posterior draw, calculate:
# - mean local decline (average across grid cells)
# - near-national decline (overall proportion of declining species)
mean_decline <- do.call(rbind, lapply(seq_along(colombia_decline), function(i) {
  draw <- colombia_decline[[i]]
  data.frame(
    draw = names(colombia_decline)[i],
    ecoregion = "Near national",
    Local = mean(draw$pointwise, na.rm = TRUE),
    Regional = draw$total,
    nsp = draw$nsp)}))


# Convert to long format for plotting
mean_long <- tidyr::pivot_longer(mean_decline, cols = c(Local, Regional),
                                 names_to = "scale", values_to = "decline")
#
# ------------------------------------------------------------------
# Compute decline summaries for each ecoregion across posterior draws
#
ecoregion_mean_decline <- rbindlist(
  lapply(names(ecoregion_decline), function(draw_name) {
    draw <- ecoregion_decline[[draw_name]]
    data.table(
      draw = draw_name,
      ecoregion = names(draw),
      Local = sapply(draw, function(x) mean(x$pointwise, na.rm = TRUE)),
      Regional = sapply(draw, function(x) x$total),
      nsp = sapply(draw, function(x) x$nsp))}))

# Clean and shorten ecoregion names for plotting
ecoregion_mean_decline$ecoregion <- gsub(" forests", "", ecoregion_mean_decline$ecoregion)
ecoregion_mean_decline$ecoregion <- 
  recode(ecoregion_mean_decline$ecoregion,
         "Apure-Villavicencio dry"         = "Villavicencio dry",
         "Cauca Valley montane"            = "Cauca montane",
         "Caqueta moist"                   ="Caquetá moist",   
         "Cordillera Oriental montane"     = "EC montane",
         "Eastern Cordillera Real montane" = "CC montane",
         "Magdalena Valley dry"            = "Magdalena dry",
         "Magdalena Valley montane"        = "Magdalena montane",
         "Northern Andean páramo"          = "Andean páramo",
         "Northwest Andean montane"        = "WC montane")

# Convert ecoregion summaries to long format for visualization
ecoregion_long <- melt(ecoregion_mean_decline,
                       id.vars = c("draw", "ecoregion", "nsp"),
                       measure.vars = c("Local", "Regional"),
                       variable.name = "scale",
                       value.name = "decline")

# Ensure consistent scale ordering
ecoregion_long$scale <- factor(ecoregion_long$scale, levels = c("Local", "Regional"))

# ------------------------------------------------------------------
# Combine and order ecoregions for plotting
decline_all <- rbind(ecoregion_long)

# Order ecoregions according to regional decline
regional_order <- decline_all %>%
  filter(scale == "Regional") %>%
  arrange(desc(decline)) %>%
  pull(ecoregion)

ecoregion_levels <- unique(c(regional_order, 
                             decline_all %>% filter(scale == "Local") %>% pull(ecoregion)))

decline_all <- decline_all %>%
  mutate(ecoregion = factor(ecoregion, levels = ecoregion_levels))

# ------------------------------------------------------------------
# Visualization of biodiversity loss distributions (Figure 3)

colors_scales <- c("Local"= "#762a83",
                   "Ecoregion" = "#1b7837",
                   "Near national" = "yellow")

# Rename scale labels for clarity in the figure
decline_all <- decline_all %>%
  rename(Scale = scale) %>%                 
  mutate(Scale = recode(Scale,             
                        "Regional" = "Ecoregion"))

# Plot posterior density distributions of biodiversity loss
fig_3 <- ggplot(decline_all, aes(x = decline, y = reorder(ecoregion, -decline), fill = Scale)) +
  geom_density_ridges(alpha = 0.7, scale = 1.5)  +
  geom_vline(xintercept = 0.86, linetype = "dashed", color = "gray30", linewidth = 0.7) +
  scale_fill_manual(values = colors_scales) +
  labs(title = NULL, x = "Proportion of declining species", y = NULL) +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(size = 11, color = "black"),
              axis.text.y = element_text(size = 11, color = "black"), 
              axis.title.x = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 13)) +
  scale_x_continuous(limits=c(0.4, 1)) +
  stat_density(geom = "line", position = "identity", aes(x = decline), kernel = "gaussian")

# This workflow summarizes posterior estimates of biodiversity loss
# across spatial scales and produces the density ridge visualization
# used in Figure 3.

ggsave("./fig_3.pdf", plot = fig_3, width = 7, height = 6, units = "in",        
       dpi = 300)
