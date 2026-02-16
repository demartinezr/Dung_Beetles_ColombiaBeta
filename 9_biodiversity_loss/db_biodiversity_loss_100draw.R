setwd("C:/Users/PC/Dropbox/CO_DBdata")
library(sf)
library(dplyr)
library(tidyr)
library(data.table)
library(purrr)
library(ggplot2)
library(ggridges)
library(patchwork)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)

################################################################################
######## get ecoregions predictions from species predictions 100 draws #########
db_predictions <- readRDS("./species_predictions_100.rds")
# prediction_data <- bind_rows(db_predictions)

# prediction_data <- prediction_data %>%
#   select(-starts_with("abun__draw_"))

# saveRDS(prediction_data, "./prediction_data_100.rds")

  # Study area ecoregions
    study_area <- readRDS("./SIG/WWF_terrestrial_ecoregions.rds")
    study_area <- st_transform(study_area, crs = 4326)
    study_area$ECO_NAME[11] <- "Caquetá moist forests"
#    study_area <- st_make_valid(study_area)
#    st_is_valid(study_area)
  
  # Function to extract abundance predictions by ecoregion
    ec_points <- vector("list", length = nrow(study_area))
    names(ec_points) <- study_area$ECO_NAME  # Use ecoregion names
  
    for (i in seq_len(nrow(study_area))) {
      cat("ecoregion:", i, "of", nrow(study_area), "\n") 
      polygon <- study_area[i, ]
    # Filter points of all species that fall within the ecoregion
    points_in_ecoregion <- do.call(rbind, lapply(db_predictions, function(df) {
      st_filter(df, polygon)
    }))
    # Store the result in the list
    ec_points[[i]] <- points_in_ecoregion
     }
# saveRDS(ec_points, "./Analysis/mean_abundance/ecoregions_predictions_100_mod.rds")
  
  # from list to sf dataframe ecoregion predictions
  ecoregions_predictions <- readRDS("./Analysis/mean_abundance/ecoregions_predictions_100_mod.rds")
    # ecoregions names in each sf dataframe
    ecoregions_predictions <- map2(
      ecoregions_predictions, 
      names(ecoregions_predictions), 
      ~ mutate(.x, ecoregions = .y)
    )
    # from list to sf dataframe
    ecoregions_all <- bind_rows(ecoregions_predictions)
#    saveRDS(ecoregions_all, "./ecoregions_100_mod.rds")

###############################################################################
#-------- Compute biodiversity loss across spatial scales and posterior
    
ecoregions_all <- as.data.frame(readRDS("./ecoregions_100_mod.rds"))
    
    paramo_sp <- c(
      "Canthidium_sp._07H", "Canthidium_sp._3", "Canthidium_sp._30H", "Canthon_arcabuquensis",
      "Canthon_sp._18H", "Canthon_sp._19H", "Canthon_sp._2", "Canthon_sp._20H",
      "Cryptocanthon_altus", "Cryptocanthon_foveatus", "Cryptocanthon_mailinae", "Deltochilum_cristinae",
      "Deltochilum_hypponum", "Deltochilum_mexicanum", "Deltochilum_sp._18H", "Deltochilum_sp._35H",
      "Dichotomius_gr._satanas", "Dichotomius_inachoides", "Dichotomius_quinquelobatus", "Eurysternus_marmoreus",
      "Homocopris_achamas", "Ontherus_brevicollis", "Ontherus_incisus", "Ontherus_lunicollis",
      "Onthophagus_curvicornis", "Onthophagus_sp._3", "Oxysternon_conspicillatum", "Uroxys_caucanus",
      "Uroxys_coarctatus", "Uroxys_cuprescens", "Uroxys_elongatus", "Uroxys_gr._pauliani",
      "Uroxys_sp._08H", "Uroxys_sp._10", "Uroxys_sp._11", "Uroxys_sp._15H",
      "Uroxys_sp._2"
    )

    ecoregions_all <- ecoregions_all %>%
      filter(
        ecoregions != "Northern Andean páramo" |
          (ecoregions == "Northern Andean páramo" & scientificName %in% paramo_sp)
      )
#    saveRDS(ecoregions_all, "./ecoregions_100_mod_P.rds")
    ecoregions_all <- readRDS("./ecoregions_100_mod_P.rds")
    # function to compute biodiversity loss  
    source("F:/repositorio/Dung_Beetles_ColombiaBeta/9_biodiversity_loss/db_compute_loss.R")
    draws <-  30 * c(1:100) # 100 draws
#    draws <- 300 * c(1:10) # 10 draws
#    draws <- 15 * c(1:200)
    unique_ecoregions <- unique(ecoregions_all$ecoregions)
  
  # Function to compute biodiversity loss across spatial scales and posterior
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
      
      eco_result$ecoregion <- eco
    ecoregion_results[[eco]] <- eco_result
    }
  
    ecoregion_ratios_list[[draw_col]] <- ecoregion_results
    rm(forest_draw, forest_draw_dt, forest_draw_wide, 
       pasture_draw, pasture_draw_dt, pasture_draw_wide, data_draw)
    gc()
  }

  saveRDS(cell_ratios_list, "./diversity_loss/cell_ratios_list_100_mod_10.rds")
  saveRDS(ecoregion_ratios_list, "./diversity_loss/ecoregion_ratios_list_100_mod_10.rds")
  saveRDS(colombia_ratios_list, "./diversity_loss/colombia_ratios_list_100_mod_10.rds")
  
################################################################################
#---------------------ecoregions sensitivity approach

# data  
  # local scale 2 x 2km
    cell_ratios_list <- readRDS("./diversity_loss/cell_ratios_list_100_mod.rds")
    cell_ratios_df <- bind_rows(cell_ratios_list, .id = "draw_col")
    cell_ratios_df <- cell_ratios_df %>%
      group_by(draw_col) %>%
      summarise(avg_logratio = mean(med_logratio, na.rm = T))
  # ecoregion scale
    ecoregion_ratios_list <- readRDS("./diversity_loss/ecoregion_ratios_list_100_mod.rds")
    ecoregion_ratios_df <- bind_rows(lapply(ecoregion_ratios_list, function(draw) {
      bind_rows(draw, .id = "ecoregion")
    }), .id = "draw_col")
    ecoregion_ratios_df <- ecoregion_ratios_df %>%
      group_by(draw_col) %>%
      summarise(avg_logratio = mean(med_logratio, na.rm = T))
  # near national scale
    colombia_ratios_list <- readRDS("./diversity_loss/colombia_ratios_list_100_mod.rds")
    colombia_ratios_df <- bind_rows(lapply(colombia_ratios_list, function(draw) {
      bind_rows(draw) %>% mutate(ecoregion = "Near national")
    }), .id = "draw_col")

###########################
## Colombia-wide comparison
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

    #relative difference
    exp(mean(colombia_ratios_df$med_logratio))/exp(mean(cell_ratios_df$med_logratio, na.rm = T))
    
    col_ratio_med <- colombia_ratios_df %>%
      select(draw_col= draw_col, col_med_logratio = med_logratio)
      
    col_diff_rel <- cell_ratios_df %>%
      group_by(draw_col) %>%
      summarise(cell_med_logratio = mean(med_logratio, na.rm = TRUE)) %>%
      left_join(col_ratio_med, by = "draw_col") %>%
      mutate(rel_dif = exp(col_med_logratio) / exp(cell_med_logratio))
    
    col_diff_summary <- col_diff_rel %>%
      summarise(avg_diff_rel = round(mean(rel_dif), 2),
                CI_5 = round(quantile(rel_dif, probs = 0.05, na.rm = TRUE), 2),
                CI_95 = round(quantile(rel_dif, probs = 0.95, na.rm = TRUE), 2))
    
# figure 2
    
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

     ggsave("./abundance_plot.jpg", plot = abundance_plot, width = 6.5, height = 4, units = "in",        
            dpi = 300, device = "jpeg")
     
########
#-- maps 
  # log-ratios map
    colombia <- ne_countries(scale = "medium", returnclass = "sf") %>%
      dplyr::filter(name == "Colombia")

    
ratio_map <- ggplot() +
      geom_sf(data = colombia, fill = "gray95", color = "black") +  
      geom_raster(data = cell_ratios_df, aes(x = lon, y = lat, fill = med_logratio)) +
      scale_fill_viridis_c(option= "B") +
      theme_bw() +
      labs(title = "a", x= NULL, y = NULL, fill = "Log ratio")  +
      theme(
        plot.title = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) + 
    scale_x_continuous(breaks = seq(-78, -68, by = 5))


  # species richness map
richness_map <- ggplot() +
  geom_sf(data = colombia, fill = "gray95", color = "black") +
  geom_raster(data = cell_ratios_df, aes(x = lon, y = lat, fill = n)) +
  scale_fill_viridis_c(option= "B") +
  theme_bw() +
  labs(title = "b", x= NULL, y = NULL, fill = "Number\n of species") +
  theme(
    plot.title = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_blank(), 
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)) + 
  scale_x_continuous(breaks = seq(-78, -68, by = 5))

log_maps <- ratio_map + richness_map

ggsave("./log_maps.jpeg", plot = log_maps, width = 8.5, height = 4, units = "in",        
       dpi = 300, device = "jpeg")
#########################################################
# ecoregions comparison across the posterior distribution

# data  
# local scale 2 x 2km
cell_ratios_eco <- bind_rows(cell_ratios_list, .id = "draw_col")
cell_ratios_eco <- cell_ratios_eco %>%
  group_by(draw_col, ecoregions) %>%
  summarise(avg_logratio = mean(med_logratio, na.rm = T))
# ecoregion scale
ecoregion_ratios_eco <- bind_rows(lapply(ecoregion_ratios_list, function(draw) {
  bind_rows(draw, .id = "ecoregion")
}), .id = "draw_col")

  # Mean ecoregion ratios
    ecoregion_mean <- ecoregion_ratios_eco %>%
      group_by(ecoregion) %>%
      summarise(med_logratio_mean = round(mean(med_logratio, na.rm = T),2),
                CI_5 = round(quantile(med_logratio, probs = 0.05, na.rm = TRUE),2),
                CI_95 = round(quantile(med_logratio, probs = 0.95, na.rm = TRUE),2))
    
  # 25th ecoregion ratios
    ecoregion_p25 <- ecoregion_ratios_eco %>%
      group_by(ecoregion) %>%
      summarise(p25_logratio_mean = round(exp(mean(p_25_logratio, na.rm = T)),2),
                CI_5 = round(quantile(p_25_logratio, probs = 0.05, na.rm = TRUE),2),
                CI_95 = round(quantile(p_25_logratio, probs = 0.95, na.rm = TRUE),2))
  
  # 75th ecoregion ratios
    ecoregion_p75 <- ecoregion_ratios_eco %>%
      group_by(ecoregion) %>%
      summarise(p_75_logratio_mean = round(mean(p_75_logratio, na.rm = T),2),
                CI_5 = round(quantile(p_75_logratio, probs = 0.05, na.rm = TRUE),2),
                CI_95 = round(quantile(p_75_logratio, probs = 0.95, na.rm = TRUE),2))
    
  # Merge ecoregion ratios with Pan-Colombia ratios
  ratios_df <- bind_rows(ecoregion_ratios_eco, colombia_ratios_df)
  
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
  
  # Region short names
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
  
    # Create a color palette function based on RdBu (red-blue)
    colors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(13))
    col_vec <- setNames(c(colors, "yellow"), regions)    
    
# Plot posterior distributions for metrics
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
    
    plot_25 / plot_mean / plot_75
    
##############################    
# Plot relative difference
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
#    relative_diff$med_ratio_dif[relative_diff$ecoregion == "Colombia"] <- NA
#    relative_diff$p_25_ratio_dif[relative_diff$ecoregion == "Colombia"] <- NA
#    relative_diff$p_75_ratio_dif[relative_diff$ecoregion == "Colombia"] <- NA
#    relative_diff <- relative_diff %>%
#      mutate(alpha_val = ifelse(ecoregion == "Colombia", 0.1, 1))
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
#------------------- Plot Relative difference cummulative
    
  # get summaries for metrics by ecoregions and near national sensitivities
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
#      filter(!ecoregion %in% c("Near national")) %>% 
      arrange(med_logratio) %>%
      mutate(ratio_exp = exp(med_logratio),
             CI_L = exp(logratio_5),
             CI_U = exp(logratio_95),
             ecoregion_order = row_number())  
    
    colombia_ratios_df
    ratios_summary_df
    
  # select draws
  draws <- paste0("abun__draw_", 30 * c(1: 100))
  
  # random sequences
  set.seed(123)
  sequences <- replicate(1000, sample(unique(ecoregion_ratios_eco$ecoregion)), simplify = FALSE)
  
  # loop for postrior distribution metrics
  
  results_all <- list()
  
  for (draw in draws) {
    cat("draw:", draw, "\n")
    
    col_draw <- colombia_ratios_df[colombia_ratios_df$draw_col == draw, ]
    ecoregion_data <- ecoregion_ratios_eco[ecoregion_ratios_eco$draw_col == draw, ]
    
    for (j in seq_along(sequences)) {
      cat("sequence", j, "\n")
      sequence <- sequences[[j]]
      
      metrics <- list(
        med_logratio = list(ref = col_draw$med_logratio, eco = ecoregion_data$med_logratio, metric = "median"),
        p_25_logratio = list(ref = col_draw$p_25_logratio, eco = ecoregion_data$p_25_logratio, metric = "p25"),
        p_75_logratio = list(ref = col_draw$p_75_logratio, eco = ecoregion_data$p_75_logratio, metric = "p75")
      )
      
      for (m in names(metrics)) {
        #sens_log <- metrics[[m]]$eco[match(sequence, ecoregion_data$ecoregion)]
        #acc_log_sens <- sapply(seq_along(sens_log), function(i) mean(sens_log[1:i], na.rm = TRUE))
        #ref_log <- mean(sens_log, na.rm = TRUE)
        #rel_diff <- exp(ref_log - acc_log_sens)
        sens <- exp(metrics[[m]]$eco[match(sequence, ecoregion_data$ecoregion)])
        acc_sens <- sapply(seq_along(sens), function(i) mean(sens[1:i]))
        ref_sens <- mean(exp(metrics[[m]]$eco), na.rm = TRUE)
        rel_diff <- ref_sens / acc_sens

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
  
  results_df <- bind_rows(results_all)
  
  # summary by metrics
  summary_df <- results_df %>%
    group_by(ecoregion, metric) %>%
    summarise(
      mean_diff = mean(relative_difference),
      p_10 = quantile(relative_difference, probs = 0.1),
      p_90 = quantile(relative_difference, probs = 0.9),
      .groups = "drop")
  
  # Plotting metrics
  med_cum <- summary_df %>% filter(metric == "median")
  p25_cum <- summary_df %>% filter(metric == "p25")
  p75_cum <- summary_df %>% filter(metric == "p75")
  
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

    fig_4 <- (plot_25 | plot_25_dif | plot_p25_cum) / 
             (plot_mean | plot_mean_dif | plot_med_cum) / 
             (plot_75 | plot_75_dif | plot_p75_cum) + 
             plot_layout(widths = c(1.5, 1, 1), heights = c(1, 1, 1))
          
  ggsave("./fig_4.pdf", plot = fig_4, width = 9, height = 9, units = "in",        
         dpi = 300)
  
