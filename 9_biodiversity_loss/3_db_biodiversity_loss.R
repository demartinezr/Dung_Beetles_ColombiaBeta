# ------------------------------------------------------------------------------
# Abundance-based sensitivity to forest–pasture conversion across spatial scales
#
# This script quantifies community sensitivity to forest conversion using
# abundance-based metrics derived from posterior predictions of the species
# abundance model. Sensitivity is calculated as the ratio of predicted abundance
# in forest relative to pasture for each species and summarized across spatial
# scales.
#
# The workflow:
# - Extracts predicted abundances within each ecoregion.
# - Computes community sensitivity metrics at local (2 × 2 km cells),
#    ecoregional, and near-national scales across posterior draws.
# 1) Summarizes these metrics to evaluate how biodiversity responses change
#    with spatial aggregation (Fig. 2b).
# 2) Examines variation in sensitivity among ecoregions and relative to the
#    near-national benchmark (Fig. 4a, d, g).
# 3) Severity of biodiversity loss at the near-national scale relative to 
#    individual ecoregions (Fig. 4b, e, h).
# 4) Assesses how estimates of biodiversity loss change as ecoregions are
#    sequentially pooled, illustrating potential underestimation when analyses
#    are restricted to few regions (Fig. 4c, f, i).
# ------------------------------------------------------------------------------

setwd("C:/Users/PC/Dropbox/CO_DBdata")
library(sf)
library(dplyr)
library(tidyr)
library(data.table)
library(purrr)
library(ggplot2)
library(ggridges)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)

# ------------------------------------------------------------------------------
# - Extract predicted abundances by ecoregion
#
# This section organizes posterior abundance predictions from the species model
# into ecoregion-level datasets. Predictions are spatially filtered using WWF
# terrestrial ecoregion polygons so that each prediction point is assigned to
# its corresponding ecoregion. The resulting datasets retain all posterior
# draws and species information, allowing subsequent sensitivity analyses to
# be computed at local, ecoregional, and near-national scales.
# ------------------------------------------------------------------------------

# Load posterior abundance predictions for all species (100 posterior draws)
  db_predictions <- readRDS("./species_predictions_100.rds")

# Study area ecoregions (WWF terrestrial ecoregions intersecting the study region)
    study_area <- readRDS("./SIG/WWF_terrestrial_ecoregions.rds")
    study_area <- st_transform(study_area, crs = 4326)

# Standardize one ecoregion name to ensure consistent matching
    study_area$ECO_NAME[11] <- "Caquetá moist forests"
    study_area <- st_make_valid(study_area)
    st_is_valid(study_area)
  
# Create a list to store predictions belonging to each ecoregion
    ec_points <- vector("list", length = nrow(study_area))
    names(ec_points) <- study_area$ECO_NAME  # Use ecoregion names

# Spatially filter prediction points so that each point is assigned
# to the ecoregion polygon in which it occurs 
    for (i in seq_len(nrow(study_area))) {
      cat("ecoregion:", i, "of", nrow(study_area), "\n") 
      polygon <- study_area[i, ]

# Extract prediction points from all species that fall inside the polygon
    points_in_ecoregion <- do.call(rbind, lapply(db_predictions, function(df) {
      st_filter(df, polygon)
    }))
    # Store the result in the list
    ec_points[[i]] <- points_in_ecoregion
    }
# Save ecoregion-specific prediction datasets (optional intermediate output)
# saveRDS(ec_points, "./Analysis/mean_abundance/ecoregions_predictions_100_mod.rds")
  
# ------------------------------------------------------------------------------
# Combine ecoregion prediction datasets
# ------------------------------------------------------------------------------
#
# Load previously saved ecoregion predictions
  ecoregions_predictions <- readRDS("./Analysis/mean_abundance/ecoregions_predictions_100_mod.rds")
#
# Add ecoregion identifiers to each dataset
    ecoregions_predictions <- map2(
      ecoregions_predictions, 
      names(ecoregions_predictions), 
      ~ mutate(.x, ecoregions = .y))
#
# Merge all ecoregion predictions into a single spatial dataframe
    ecoregions_all <- bind_rows(ecoregions_predictions)
#   saveRDS(ecoregions_all, "./ecoregions_100_mod_P.rds")
#
# ------------------------------------------------------------------------------
# - Compute abundance-based sensitivity to forest conversion across spatial scales
#
# Using posterior abundance predictions organized by ecoregion, this section
# calculates community sensitivity to forest–pasture conversion. Sensitivity is
# defined as the ratio of predicted abundance in forest relative to pasture for
# each species and summarized across communities.
#
# For each posterior draw, sensitivity metrics are computed at three spatial
# scales:
# 1) Local scale (2 × 2 km grid cells)
# 2) Ecoregional scale (aggregating cells within each ecoregion)
# 3) Near-national scale (all cells combined)
#
# The resulting metrics are stored for each posterior draw and later used to
# quantify biodiversity responses to forest conversion (Fig. 2b) and spatial
# variation in sensitivity across ecoregions (Fig. 3).
  

# Load ecoregion-level prediction dataset
    ecoregions_all <- readRDS("./ecoregions_100_mod_P.rds")
    
# Load functions used to compute community sensitivity metrics  
    source("F:/repositorio/Dung_Beetles_ColombiaBeta/9_biodiversity_loss/db_compute_loss.R")
    
# Select posterior draws used in the analysis
    draws <-  30 * c(1:100)

# Identify unique ecoregions in the dataset
    unique_ecoregions <- unique(ecoregions_all$ecoregions)
  
# Initialize lists to store sensitivity metrics for each posterior draw
  cell_ratios_list <- list()
  colombia_ratios_list <- list()
  ecoregion_ratios_list <- list()

  for (i in seq_along(draws)) {
      cat("draw:", i, "of", length(draws), "\n")
    # Extract data corresponding to the current posterior draw
      draw_col <- paste0("abun__draw_", draws[i])
    # Filter data and add ID cell
    data_draw <- ecoregions_all %>%
      select(ecoregions, lon, lat, scientificName, pasture, all_of(draw_col)) %>%
      mutate(cell_id = row_number()) %>%
      select(cell_id, everything())
    row.names(data_draw) <- NULL
    # Split predictions into forest and pasture communitie
    forest_draw <- data_draw %>% filter(pasture == 0) %>% select(-pasture)
    pasture_draw <- data_draw %>% filter(pasture == 1) %>% select(-pasture)
    # Convert to community matrices (sites × species)
    forest_draw_dt <- as.data.table(forest_draw)
    forest_draw_wide <- dcast(forest_draw_dt, ecoregions + lon + lat ~ scientificName, 
                              value.var = draw_col, fill = 0)
    
    pasture_draw_dt <- as.data.table(pasture_draw)
    pasture_draw_wide <- dcast(pasture_draw_dt, ecoregions + lon + lat ~ scientificName, 
                               value.var = draw_col, fill = 0)
    # Compute sensitivity metrics at local and near-national scales
    cell_ratios_list[[draw_col]] <- get_avg_cell_ratios(forest_draw_wide, pasture_draw_wide, 
                                                        cutoff_type = "absolute", cutoff = 1)
    colombia_ratios_list[[draw_col]] <- get_regional_ratios(forest_draw_wide, pasture_draw_wide, 
                                                            cutoff_type = "absolute", cutoff = 1, 
                                                            cell_positions = NULL)
    # Compute sensitivity metrics for each ecoregion
    ecoregion_results <- list()
    
    for (j in seq_along(unique_ecoregions)) {
      cat(j, "\n") 
      
      eco <- unique_ecoregions[j]
      cell_positions <- which(forest_draw_wide$ecoregions == eco)
      eco_result <- get_regional_ratios(forest_draw_wide, pasture_draw_wide, 
                                        cutoff_type = "absolute", cutoff = 1, 
                                        cell_positions = cell_positions)
      
      eco_result$ecoregion <- eco
    ecoregion_results[[eco]] <- eco_result
    }
  
    ecoregion_ratios_list[[draw_col]] <- ecoregion_results
    rm(forest_draw, forest_draw_dt, forest_draw_wide, 
       pasture_draw, pasture_draw_dt, pasture_draw_wide, data_draw)
    gc()
  }
  
  # Save posterior sensitivity metrics for downstream analyses
  saveRDS(cell_ratios_list, "./diversity_loss/cell_ratios_list_100_mod_10.rds")
  saveRDS(ecoregion_ratios_list, "./diversity_loss/ecoregion_ratios_list_100_mod_10.rds")
  saveRDS(colombia_ratios_list, "./diversity_loss/colombia_ratios_list_100_mod_10.rds")
# 
# ------------------------------------------------------------------------------
#  1) Median community-level sensitivity to forest-pasture conversion
#
# This section evaluates how biodiversity responses to forest–pasture
# differences change when results are aggregated across spatial scales.
# The analysis summarizes median forest–pasture log-ratios across posterior
# draws at three spatial levels: local grid cells, ecoregions, and near-national
# scale (Fig. 2b).

# data  
  # Local scale (2 × 2 km grid cells)
    # Median forest–pasture log-ratios per grid cell for each posterior draw
    cell_ratios_list <- readRDS("./diversity_loss/cell_ratios_list_100_mod.rds")
    cell_ratios_df <- bind_rows(cell_ratios_list, .id = "draw_col")
    
    # Average median log-ratio across cells for each posterior draw
    cell_ratios_df <- cell_ratios_df %>%
      group_by(draw_col) %>%
      summarise(avg_logratio = mean(med_logratio, na.rm = T))
    
  # Ecoregion scale
  # Median log-ratios aggregated within ecoregions
    ecoregion_ratios_list <- readRDS("./diversity_loss/ecoregion_ratios_list_100_mod.rds")
    ecoregion_ratios_df <- bind_rows(lapply(ecoregion_ratios_list, function(draw) {
      bind_rows(draw, .id = "ecoregion")
    }), .id = "draw_col")
    
    # Average median log-ratio across ecoregions for each posterior draw
    ecoregion_ratios_df <- ecoregion_ratios_df %>%
      group_by(draw_col) %>%
      summarise(avg_logratio = mean(med_logratio, na.rm = T))
    
  # Near-national scale
  # Log-ratios aggregated across all ecoregions
    colombia_ratios_list <- readRDS("./diversity_loss/colombia_ratios_list_100_mod.rds")
    colombia_ratios_df <- bind_rows(lapply(colombia_ratios_list, function(draw) {
      bind_rows(draw) %>% mutate(ecoregion = "Near national")
    }), .id = "draw_col")

# ------------------------------------------------------------------------------
# Colombia-wide comparison
    
    # Mean forest–pasture log-ratios and abundance ratios across posterior draws
    
  # local scale 2x2km metrics
    mean(cell_ratios_df$med_logratio, na.rm = T)
    mean(exp(cell_ratios_df$med_logratio), na.rm = T)
    exp(mean(cell_ratios_df$med_logratio, na.rm = T))
    # ecoregions metrics
    mean(ecoregion_ratios_df$med_logratio, na.rm = T)
    mean(exp(ecoregion_ratios_df$med_logratio))
    exp(mean(ecoregion_ratios_df$med_logratio))
  # near national metrics  
    mean(colombia_ratios_df$med_logratio, na.rm = T)
    mean(exp(colombia_ratios_df$med_logratio))
    exp(mean(colombia_ratios_df$med_logratio))

  # Credible intervals of abundance ratios
  # IC local scale
    quantile(cell_ratios_df$med_logratio, probs = c(0.05, 0.95), na.rm = TRUE)
    round(quantile(exp(cell_ratios_df$med_logratio), probs = c(0.05, 0.95), na.rm = TRUE),2)
    round((exp(quantile(cell_ratios_df$med_logratio, probs = c(0.05, 0.95), na.rm = TRUE)) - 1) * 100, 1)
  # IC ecoregions
    quantile(ecoregion_ratios_df$med_logratio, probs = c(0.05, 0.95))
    round(quantile(exp(ecoregion_ratios_df$med_logratio), probs = c(0.05, 0.95)), 2) # absolute
    round((exp(quantile(ecoregion_ratios_df$med_logratio, probs = c(0.05, 0.95), na.rm = TRUE)) - 1) * 100, 1) # percent 
  # IC near national
    quantile(colombia_ratios_df$med_logratio, probs = c(0.05, 0.95))
    round(quantile(exp(colombia_ratios_df$med_logratio), probs = c(0.05, 0.95)), 2) # absolute
    round((exp(quantile(colombia_ratios_df$med_logratio, probs = c(0.05, 0.95), na.rm = TRUE)) - 1) * 100, 1) # percent 

# Relative difference between near-national and local-scale abundance ratios
    exp(mean(colombia_ratios_df$med_logratio))/exp(mean(cell_ratios_df$med_logratio, na.rm = T))
    
    col_ratio_med <- colombia_ratios_df %>%
      select(draw_col= draw_col, col_med_logratio = med_logratio)

# Compute relative difference per posterior draw          
    col_diff_rel <- cell_ratios_df %>%
      group_by(draw_col) %>%
      summarise(cell_med_logratio = mean(med_logratio, na.rm = TRUE)) %>%
      left_join(col_ratio_med, by = "draw_col") %>%
      mutate(rel_dif = exp(col_med_logratio) / exp(cell_med_logratio))

# Summary statistics of relative differences    
    col_diff_summary <- col_diff_rel %>%
      summarise(avg_diff_rel = round(mean(rel_dif), 2),
                CI_5 = round(quantile(rel_dif, probs = 0.05, na.rm = TRUE), 2),
                CI_95 = round(quantile(rel_dif, probs = 0.95, na.rm = TRUE), 2))
    
# Figure 2b
# Distribution of forest–pasture log-ratios across spatial scales
    
abundance_scales <- bind_rows(
      cell_ratios_df %>%
        mutate(scale = "Local"),
      
      ecoregion_ratios_df %>%
        mutate(scale = "Ecoregion"),
      
      colombia_ratios_df %>%
        mutate(scale = "Near national"))
    
abundance_scales$scale <- factor(abundance_scales$scale,
                                    levels = c("Local", "Ecoregion", "Near national"))
    
colors_scales <- c("Local"= "#762a83",
                       "Ecoregion" = "#1b7837",
                       "Near national" = "yellow")
        
fig_2b <- ggplot(abundance_scales, aes(x = reorder(scale, avg_logratio), y = avg_logratio, fill = scale)) +
                     geom_jitter(width = 0.2, size = 1, alpha = 0.4, color="gray80") +
                     geom_violin(trim = FALSE, alpha = 0.5) +
                     geom_boxplot(width = 0.1, outlier.shape = NA) +
                     scale_fill_manual(values = colors_scales) +
                     labs(title = "b", x = NULL, y = "Sensitivity", fill = "Scale") +
                     theme_bw() +
                     theme(legend.position = "none",
                           #legend.position.inside  = c(0.99, 0.02),  
                           #legend.justification = c("right", "bottom"),
                           legend.direction = "horizontal",       
                           #legend.box.background = element_rect(color = "black", fill = "white", linewidth = 0.5),
                           #legend.box.margin = margin(2, 2, 2, 2),
                           plot.title = element_text(size = 12, face = "bold"),
                           legend.text = element_text(size = 12),
                           legend.title = element_text(size = 13),
                           axis.title.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.y = element_text(size = 13, color = "black"),
                           axis.text.y = element_text(size = 12, color = "black"),
                           axis.text.x = element_text(size = 12, color = "black"))

################################################################################
# 2) Community sensitivity to forest–pasture conversion (change in abundance) 
# for each ecoregion in absolute terms
#
# This section evaluates how community sensitivity to forest–pasture conversion
# varies among ecoregions. Posterior log-ratios of abundance (forest/pasture)
# are summarized for each ecoregion and compared with the near-national
# distribution. Sensitivity is examined across the species response spectrum,
# including low (25th percentile), median (50th percentile), and high
# sensitivity species (75th percentile), forming the basis of Fig.4a, d, g
# ------------------------------------------------------------------------------

# data

# Local scale (2 × 2 km grid cells)
# Average median log-ratios per ecoregion for each posterior draw
cell_ratios_eco <- bind_rows(cell_ratios_list, .id = "draw_col")
cell_ratios_eco <- cell_ratios_eco %>%
  group_by(draw_col, ecoregions) %>%
  summarise(avg_logratio = mean(med_logratio, na.rm = T))

# Ecoregion scale
# Posterior log-ratio summaries computed directly at the ecoregion level
ecoregion_ratios_eco <- bind_rows(lapply(ecoregion_ratios_list, function(draw) {
  bind_rows(draw, .id = "ecoregion")
}), .id = "draw_col")

# Summary statistics of ecoregion sensitivity

    # Mean (median sensitivity of species)
    ecoregion_mean <- ecoregion_ratios_eco %>%
      group_by(ecoregion) %>%
      summarise(med_logratio_mean = round(mean(med_logratio, na.rm = T),2),
                CI_5 = round(quantile(med_logratio, probs = 0.05, na.rm = TRUE),2),
                CI_95 = round(quantile(med_logratio, probs = 0.95, na.rm = TRUE),2))
    
    # Low-sensitivity species (25th percentile of species responses)
    ecoregion_p25 <- ecoregion_ratios_eco %>%
      group_by(ecoregion) %>%
      summarise(p25_logratio_mean = round(exp(mean(p_25_logratio, na.rm = T)),2),
                CI_5 = round(quantile(p_25_logratio, probs = 0.05, na.rm = TRUE),2),
                CI_95 = round(quantile(p_25_logratio, probs = 0.95, na.rm = TRUE),2))
  
    # High-sensitivity species (75th percentile of species responses)
    ecoregion_p75 <- ecoregion_ratios_eco %>%
      group_by(ecoregion) %>%
      summarise(p_75_logratio_mean = round(mean(p_75_logratio, na.rm = T),2),
                CI_5 = round(quantile(p_75_logratio, probs = 0.05, na.rm = TRUE),2),
                CI_95 = round(quantile(p_75_logratio, probs = 0.95, na.rm = TRUE),2))
    
  # Combine ecoregion and near-national posterior distributions
  ratios_df <- bind_rows(ecoregion_ratios_eco, colombia_ratios_df)
 
  # Define ecoregion order used in the figures
  regions <- c("Eastern Cordillera real montane forests",
               "Magdalena-Urabá moist forests",
               "Apure-Villavicencio dry forests",
               "Napo moist forests",
               "Cordillera Oriental montane forests",
               "Llanos",
               "Santa Marta montane forests",
               "Magdalena Valley montane forests",
               "Northern Andean páramo",
               "Caquetá moist forests",
               "Magdalena Valley dry forests",
               "Cauca Valley montane forests",
               "Northwestern Andean montane forests",
               "Near national")
  
  ratios_df$ecoregion <- factor(ratios_df$ecoregion, levels = regions)
  
  # Short labels used in the figure
  reg_short <- c("CC montane",
                 "Magdalena-Urabá moist",
                 "Villavicencio dry",
                 "Napo moist",
                 "EC montane",
                 "Llanos",
                 "Santa Marta montane",
                 "Magdalena montane",
                 "Andean páramo",
                 "Caquetá moist",
                 "Magdalena dry",
                 "Cauca montane",
                 "WC montane",
                 "Near national")
  label_vec <- setNames(reg_short, regions)
  
  # Color palette used to distinguish ecoregions in density plots
    colors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(13))
    col_vec <- setNames(c(colors, "yellow"), regions)    
    
# Plot posterior distributions of community sensitivity
# across ecoregions (Fig. 4a, d, g)
    
    
# Fig 4a. Low sensitivity species (25th percentile) 
    plot_25 <- ggplot(ratios_df, aes(x = p_25_logratio, y = reorder(ecoregion, -p_25_logratio), fill = ecoregion)) +
      geom_density_ridges(alpha = 0.7, scale = 1.5)  +
      scale_fill_manual(values = col_vec, labels = label_vec) +
      scale_y_discrete(labels = label_vec) + 
      labs(title = "a", x = "sensitivity\n(N forest/N pasture)", y = "Low sensitivity species\n (25th percentile)") +
      theme_bw() +
      theme(legend.position = "none",
        axis.text.x = element_text(size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 12, face = "bold")) +
      scale_x_continuous(limits=c(-1, 8)) +
      stat_density(geom = "line",position = "identity", aes(x = p_25_logratio), kernel = "gaussian")

# Fig 4d. Median sensitivity species (50th percentile)    
    plot_mean <- ggplot(ratios_df, aes(x = med_logratio, y = reorder(ecoregion, -med_logratio), fill = ecoregion)) +
      geom_density_ridges(alpha = 0.7, scale = 1.5)  +
      scale_fill_manual(values = col_vec, labels = label_vec) +
      scale_y_discrete(labels = label_vec) + 
      labs(title = "d", x = "sensitivity\n(N forest/N pasture)", y = "Medium sensitivity species\n (50th percentile)") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 11, color = "black"),
            axis.text.y = element_text(size = 11, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 12),
            plot.title = element_text(size = 12, face = "bold")) +
      scale_x_continuous(limits=c(-1, 8)) +
      stat_density(geom = "line", position = "identity", aes(x = med_logratio), kernel = "gaussian")

# Fig 4g. High sensitivity species (75th percentile) 
    plot_75 <- ggplot(ratios_df, aes(x = p_75_logratio, y = reorder(ecoregion, -p_75_logratio), fill = ecoregion)) +
      geom_density_ridges(alpha = 0.7, scale = 1.5)  +
      scale_fill_manual(values = col_vec, labels = label_vec) +
      scale_y_discrete(labels = label_vec) + 
      labs(title = "g", x = "sensitivity (N forest/N pasture)", y = "High sensitivity species\n (75th percentile)") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 11, color = "black"),
            axis.text.y = element_text(size = 11, color = "black"),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            plot.title = element_text(size = 12, face = "bold")) +
      scale_x_continuous(limits=c(-1, 8)) +
      stat_density(geom = "line", position = "identity", aes(x = p_75_logratio), kernel = "gaussian")

# ------------------------------------------------------------------------------
#  3) Severity of biodiversity loss at the near-national scale relative to individual ecoregions.
#
# This section quantifies how biodiversity responses in each ecoregion differ
# from the near-national estimate. For each posterior draw, ecoregion-level
# richness ratios (forest / pasture) are compared with the corresponding
# near-national ratio to obtain a relative difference metric.
#
# These relative differences are summarized across posterior draws to obtain
# mean estimates and uncertainty intervals for each ecoregion. The resulting
# distributions are then visualized using ridge density plots, showing how
# the lower (25th percentile), median, and upper (75th percentile) portions
# of the species response distribution deviate from the near-national baseline
# (Fig. 4b, e, h).
# ------------------------------------------------------------------------------
    
# Summarize relative differences across posterior draws for each ecoregion
# to obtain mean estimates and uncertainty intervals    
    relative_diff <- ratios_df %>%
      left_join(colombia_ratios_df, by = "draw_col", suffix = c("", "_col")) %>% 
      mutate(
        avg_ratio_dif = avg_ratio_col / avg_ratio,
        avg_logratio_dif  = exp(avg_logratio_col - avg_logratio),
        med_ratio_dif  = exp(med_logratio_col - med_logratio),
        p_25_ratio_dif  = exp(p_25_logratio_col - p_25_logratio),
        p_75_ratio_dif  = exp(p_75_logratio_col - p_75_logratio)) %>%
      select(ecoregion, draw_col, avg_ratio_dif, avg_logratio_dif, med_ratio_dif, p_25_ratio_dif, p_75_ratio_dif)
    
    
    relative_diff_mean <-relative_diff %>%
      group_by(ecoregion) %>%
      summarise(
        avg_ratio_dif_mean = mean(avg_ratio_dif, na.rm = TRUE),
        avg_logratio_dif_mean = mean(avg_logratio_dif, na.rm = TRUE),
        med_logratio_dif_mean = mean(med_ratio_dif, na.rm = TRUE),
        med_p5 = quantile(med_ratio_dif, probs = 0.05, na.rm = TRUE),
        med_p95 = quantile(med_ratio_dif, probs = 0.95, na.rm = TRUE),
        p_25_logratio_dif = mean(p_25_ratio_dif, na.rm = TRUE),
        p25_p5 = quantile(p_25_ratio_dif, probs = 0.05, na.rm = TRUE),
        p25_p95 = quantile(p_25_ratio_dif, probs = 0.95, na.rm = TRUE),
        p_75_logratio_dif = mean(p_75_ratio_dif, na.rm = TRUE),
        p75_p5 = quantile(p_75_ratio_dif, probs = 0.05, na.rm = TRUE),
        p75_p95 = quantile(p_75_ratio_dif, probs = 0.95, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(med_logratio_dif_mean) %>%
      mutate(ecoregion_order = row_number())
    
# Summarize the original ecoregion ratios (used to define plotting order)    
    ratios_df_mean <-ratios_df %>%
      group_by(ecoregion) %>%
      summarise(
        avg_ratio_mean = mean(avg_ratio, na.rm = TRUE),
        avg_logratio_mean = mean(avg_logratio, na.rm = TRUE),
        med_logratio_mean = mean(med_logratio, na.rm = TRUE),
        med_p5 = quantile(med_logratio, probs = 0.05, na.rm = TRUE),
        med_p95 = quantile(med_logratio, probs = 0.95, na.rm = TRUE),
        p_25_logratio = mean(p_25_logratio, na.rm = TRUE),
        p25_p5 = quantile(p_25_logratio, probs = 0.05, na.rm = TRUE),
        p25_p95 = quantile(p_25_logratio, probs = 0.95, na.rm = TRUE),
        p_75_logratio = mean(p_75_logratio, na.rm = TRUE),
        p75_p5 = quantile(p_75_logratio, probs = 0.05, na.rm = TRUE),
        p75_p95 = quantile(p_75_logratio, probs = 0.95, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(med_logratio_mean) %>%
      mutate(ecoregion_order = row_number())

# Remove the near-national category from density plots so that
# only deviations of individual ecoregions are visualized
      relative_diff_dum <- relative_diff %>%
      mutate(across(
        c(avg_ratio_dif, avg_logratio_dif, med_ratio_dif, 
          p_25_ratio_dif, p_75_ratio_dif),
        ~ ifelse(ecoregion == "Near national", NA, .)
      ))
        
    relative_diff_dum$ecoregion <- factor(
      relative_diff_dum$ecoregion,
      levels = ratios_df_mean %>%
        arrange(desc(p_25_logratio)) %>%
        pull(ecoregion))
    
# Ridge density plots showing how the 25th percentile, median,
# and 75th percentile of the response distribution differ
# from the near-national baseline across ecoregions
    
    # Fig 4b. Lower-sensitivity species (25th percentile) relative to near-national baseline
    plot_25_dif <- ggplot(relative_diff_dum, aes(x = p_25_ratio_dif, y = ecoregion, fill = ecoregion)) +
      geom_density_ridges(alpha = 0.7, scale = 1.5)  +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
      scale_fill_manual(values = col_vec, labels = label_vec) +
      scale_y_discrete(labels = label_vec) + 
      labs(title = "b", x = "sensitivity\n(N forest/N pasture)", y = "Ecoregion") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 11, color = "black"),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 12, face = "bold")) +
      scale_x_continuous(limits = c(-1, 12), breaks = c(1, 4, 8, 12, 16, 20)) +
      stat_density(geom = "line",position = "identity", aes(x = p_25_ratio_dif), kernel = "gaussian")
    
    relative_diff_dum$ecoregion <- factor(
      relative_diff_dum$ecoregion,
      levels = ratios_df_mean %>%
        arrange(desc(med_logratio_mean)) %>%
        pull(ecoregion))  
    
    # Fig. 4e. Medium-sensitivity species (median / 50th percentile) relative to near-national baseline
    plot_mean_dif <- ggplot(relative_diff_dum, aes(x = med_ratio_dif, y = ecoregion, fill = ecoregion)) +
      geom_density_ridges(alpha = 0.7, scale = 1.5)  +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
      scale_fill_manual(values = col_vec, labels = label_vec) +
      scale_y_discrete(labels = label_vec) + 
      labs(title = "e", x = "sensitivity\n(N forest/N pasture)", y = "Ecoregion") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 11, color = "black"),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 12, face = "bold")) +
      scale_x_continuous(limits = c(-1, 12), breaks = c(1, 4, 8, 12, 16, 20))+
      stat_density(geom = "line", position = "identity", aes(x = med_ratio_dif), kernel = "gaussian")
    
    relative_diff_dum$ecoregion <- factor(
      relative_diff_dum$ecoregion,
      levels = ratios_df_mean %>%
        arrange(desc(p_75_logratio)) %>%
        pull(ecoregion)) 
    
    # Fig 4g. High-sensitivity species (75th percentile) relative to near-national baseline
    plot_75_dif <- ggplot(relative_diff_dum, aes(x = p_75_ratio_dif, y = ecoregion, fill = ecoregion)) +
      geom_density_ridges(alpha = 0.7, scale = 1.5)  +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
      scale_fill_manual(values = col_vec, labels = label_vec) +
      scale_y_discrete(labels = label_vec) + 
      labs(title = "h", x = "Relative difference", y = "Ecoregion") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 11, color = "black"),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 12, face = "bold")) +
      scale_x_continuous(limits = c(-1, 12), breaks = c(1, 4, 8, 12, 16, 20)) +
      stat_density(geom = "line", position = "identity", aes(x = p_75_ratio_dif), kernel = "gaussian")

################################################################################
#  4) Regional pooling to assess underestimation by local-scale reliance
#
# This section evaluates how biodiversity loss estimates change as data from 
# ecoregions are sequentially pooled, relative to the near-national benchmark.
# It quantifies underestimation of biodiversity loss when analyses are restricted
# to a few ecoregions, capturing not only median responses but also species-level 
# distributions (25th and 75th percentiles). These cumulative plots are shown
# in Fig. 4c, f, i.
# ------------------------------------------------------------------------------

# Compute summary metrics for each ecoregion across posterior draws
# including mean, 10th and 90th percentiles, median, 25th and 75th percentiles
    ratios_summary_df <- ecoregion_ratios_eco %>%
      group_by(ecoregion) %>%
      summarise(
        avg_ratio = mean(avg_ratio, na.rm = TRUE),
        logratio_5 = quantile(med_logratio, probs = 0.1, na.rm = TRUE),
        logratio_95 = quantile(med_logratio, probs = 0.9, na.rm = TRUE),
        avg_logratio = mean(avg_logratio, na.rm = TRUE),
        med_logratio = mean(med_logratio, na.rm = TRUE),
        p_25_logratio = mean(p_25_logratio, na.rm = TRUE),
        p_75_logratio = mean(p_75_logratio, na.rm = TRUE),
        .groups = "drop") %>%
      arrange(med_logratio) %>%
      mutate(ratio_exp = exp(med_logratio),
             CI_L = exp(logratio_5),
             CI_U = exp(logratio_95),
             ecoregion_order = row_number())  

   # Reference to Colombia-wide ratios (near-national baseline)    
    colombia_ratios_df
    ratios_summary_df
    
  # Define posterior draw names
  draws <- paste0("abun__draw_", 30 * c(1: 100))
  
  # Generate 1000 random sequences of ecoregions for cumulative pooling
  set.seed(123)
  sequences <- replicate(1000, sample(unique(ecoregion_ratios_eco$ecoregion)), simplify = FALSE)
  
  # Initialize list to store cumulative result
  results_all <- list()
  # Loop through posterior draws
  for (draw in draws) {
    cat("draw:", draw, "\n")
    # Select data for current draw
    col_draw <- colombia_ratios_df[colombia_ratios_df$draw_col == draw, ]
    ecoregion_data <- ecoregion_ratios_eco[ecoregion_ratios_eco$draw_col == draw, ]
    # Loop through each random ecoregion sequence
    for (j in seq_along(sequences)) {
      cat("sequence", j, "\n")
      sequence <- sequences[[j]]
      # Define metrics for median, 25th, and 75th percentile
      metrics <- list(
        med_logratio = list(ref = col_draw$med_logratio, eco = ecoregion_data$med_logratio, metric = "median"),
        p_25_logratio = list(ref = col_draw$p_25_logratio, eco = ecoregion_data$p_25_logratio, metric = "p25"),
        p_75_logratio = list(ref = col_draw$p_75_logratio, eco = ecoregion_data$p_75_logratio, metric = "p75")
      )
      # Loop through metrics and compute cumulative relative differences
      for (m in names(metrics)) {
        # Sensitivity values for the sequence
        sens <- exp(metrics[[m]]$eco[match(sequence, ecoregion_data$ecoregion)])
        # Cumulative mean sensitivity as ecoregions are pooled
        acc_sens <- sapply(seq_along(sens), function(i) mean(sens[1:i]))
        # Reference near-national mean sensitivity
        ref_sens <- mean(exp(metrics[[m]]$eco), na.rm = TRUE)
        # Relative difference (underestimation factor)
        rel_diff <- ref_sens / acc_sens
        # Store results
        results_all[[length(results_all) + 1]] <- data.frame(
          draw = draw,
          sequence = j,
          ecoregion = seq_along(sequence),
          relative_difference = rel_diff,
          metric = metrics[[m]]$metric
        )
      }
    }
  }
 
  # Combine all results into a single dataframe
  results_df <- bind_rows(results_all)
  
  # Summarize relative differences by ecoregion and metric
  summary_df <- results_df %>%
    group_by(ecoregion, metric) %>%
    summarise(
      mean_diff = mean(relative_difference),
      p_10 = quantile(relative_difference, probs = 0.1),
      p_90 = quantile(relative_difference, probs = 0.9),
      .groups = "drop")
  
  # Split summary by metric for plotting
  med_cum <- summary_df %>% filter(metric == "median")
  p25_cum <- summary_df %>% filter(metric == "p25")
  p75_cum <- summary_df %>% filter(metric == "p75")
  
  # Plot cumulative relative differences for low, medium, and high sensitivity species
  # Black points and error bars = posterior means and 90% credible intervals
  # Grey ribbon = range around near-national benchmark (dashed line at 1)
  
  # Fig 4c. 25th percentile cumulative relative difference (low-sensitivity species)
  plot_p25_cum <- ggplot(p25_cum, aes(x = ecoregion, y = mean_diff)) +
    geom_point(size = 2.5, color = "black") +  # Puntos para la media
    geom_errorbar(aes(ymin = p_10, ymax = p_90), width = 0, color = "black", linewidth = 0.6) +  # Barras de error
    geom_ribbon(aes(ymin = 0.99, ymax = 1.05), fill = "gray", alpha = 0.5) +  # Banda gris
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Línea de referencia en 1
    scale_x_continuous(breaks = seq(1, 13, by = 2)) +  # Escala del eje X
    scale_y_continuous(limits=c(0.3, 5), breaks = seq(1, 5, by = 1)) +
    labs(title = "c", x = NULL, y = "Relative difference") +
    theme_bw() +
    theme(plot.title = element_text(size = 12, face = "bold"), 
          axis.title.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 11, color = "black"),
          axis.text.y = element_text(size = 11, color = "black"))

  # Fig 4f. 50th percentile (median) cumulative relative difference    
  plot_med_cum <- ggplot(med_cum, aes(x = ecoregion, y = mean_diff)) +
    geom_point(size = 2.5, color = "black") +  # Puntos para la media
    geom_errorbar(aes(ymin = p_10, ymax = p_90), width = 0, color = "black", linewidth = 0.6) +  # Barras de error
    geom_ribbon(aes(ymin = 0.99, ymax = 1.05), fill = "gray", alpha = 0.5) +  # Banda gris
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Línea de referencia en 1
    scale_x_continuous(breaks = seq(1, 13, by = 2)) +  # Escala del eje X
    scale_y_continuous(limits=c(0.3, 5), breaks = seq(1, 5, by = 1)) +
    labs(title = "f", x = NULL, y = "Relative difference") +
    theme_bw() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.title.x = element_text(size  = 16, face = "italic"),
          axis.text.x = element_text(size = 11, color = "black"),
          axis.text.y = element_text(size = 11, color = "black"), 
          axis.title.y = element_text(size = 11))

  # Fig 4i. 75th percentile cumulative relative difference (high-sensitivity species)    
  plot_p75_cum <- ggplot(p75_cum, aes(x = ecoregion, y = mean_diff)) +
    geom_point(size = 2.5, color = "black") + 
    geom_errorbar(aes(ymin = p_10, ymax = p_90),  width = 0, color = "black", linewidth = 0.6) +  # Barras de error
    geom_ribbon(aes(ymin = 0.99, ymax = 1.05), fill = "gray", alpha = 0.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") + 
    scale_x_continuous(breaks = seq(1, 13, by = 2)) +
    scale_y_continuous(limits=c(0.3, 5), breaks = seq(1, 5, by = 1)) +
    labs(title ="i", x = expression(N[regions]), y = "Relative difference") +
    theme_bw() +
    theme(plot.title = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 11, color = "black"), 
          axis.title.y = element_text(size  = 12))

# ------------------------------------------------------------------------------
# Assemble Figure 4: Multi-scale sensitivity of dung beetle communities
#
# This layout combines all panels showing species-level sensitivity to forest–pasture
# conversion across ecoregions and scales:
# - Left column: Low (25th percentile), medium (50th percentile), and high (75th percentile) sensitivity distributions
#   at the ecoregion scale (panels a, d, g).
# - Middle column: Relative difference compared to near-national baseline (panels b, e, h; dashed vertical line = 1).
# - Right column: Cumulative relative difference as ecoregions are sequentially pooled (panels c, f, i).
# Panel widths are adjusted to emphasize the distribution plots, and heights are uniform across rows.
  
    fig_4 <- (plot_25 | plot_25_dif | plot_p25_cum) / 
             (plot_mean | plot_mean_dif | plot_med_cum) / 
             (plot_75 | plot_75_dif | plot_p75_cum) + 
             plot_layout(widths = c(1.5, 1, 1), heights = c(1, 1, 1))
          
  ggsave("./fig_4.pdf", plot = fig_4, width = 9, height = 9, units = "in",        
         dpi = 300)
  
