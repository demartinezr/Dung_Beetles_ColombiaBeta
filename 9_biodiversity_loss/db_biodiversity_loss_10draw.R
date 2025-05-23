setwd("C:/Users/PC/Dropbox/CO_DBdata")
library(sf)
library(dplyr)
library(tidyr)
library(data.table)
library(purrr)
library(ggplot2)
library(ggridges)
library(gridExtra)

################################################################################
######## get ecoregions predictions from species predictions 100 draws #########
db_predictions <- readRDS("./species_predictions_100.rds")
# prediction_data <- bind_rows(db_predictions)

# prediction_data <- prediction_data %>%
#   select(-starts_with("abun__draw_"))

# saveRDS(prediction_data, "./prediction_data_100.rds")

  # Study area ecoregions
    study_area <- st_read("F:/Capas/America/ecoregions/ecoreg.shp")
    study_area <- st_make_valid(study_area)
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
# saveRDS(ec_points, "./Analysis/mean_abundance/ecoregions_predictions_100.rds")
  
  # from list to sf dataframe ecoregion predictions
  ecoregions_predictions <- readRDS("./Analysis/mean_abundance/ecoregions_predictions_100.rds")
    # ecoregions names in each sf dataframe
    ecoregions_predictions <- map2(
      ecoregions_predictions, 
      names(ecoregions_predictions), 
      ~ mutate(.x, ecoregions = .y)
    )
    # from list to sf dataframe
    ecoregions_all <- bind_rows(ecoregions_predictions)
    saveRDS(ecoregions_all, "./ecoregions_100.rds")

################################################################################
####### Compute biodiversity loss across spatial scales and posterior ##########
ecoregions_all <- as.data.frame(readRDS("./ecoregions_100.rds"))

  # function to compute biodiversity loss  
    source("F:/repositorio/Dung_Beetles_ColombiaBeta/9_biodiversity_loss/db_compute_loss.R")
    draws <-  30 * c(1:100) # 100 draws
#    draws <- 300 * c(1:10) # 10 draws
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
      
      eco_result$ecoregion <- eco  # Agregar la ecorregión a los resultados
    ecoregion_results[[eco]] <- eco_result
    }
  
    ecoregion_ratios_list[[draw_col]] <- ecoregion_results
    rm(forest_draw, forest_draw_dt, forest_draw_wide, 
       pasture_draw, pasture_draw_dt, pasture_draw_wide, data_draw)
    gc()
  }

  saveRDS(cell_ratios_list, "./diversity_loss/cell_ratios_list_100.rds")
  saveRDS(ecoregion_ratios_list, "./diversity_loss/ecoregion_ratios_list_100.rds")
  saveRDS(colombia_ratios_list, "./diversity_loss/colombia_ratios_list_100.rds")
  
################################################################################
#################### ecoregions sensitivity approach ###########################

# load data  
  # local scale 2km
    cell_ratios_list <- readRDS("./diversity_loss/cell_ratios_list_100.rds")
    cell_ratios_df <- bind_rows(cell_ratios_list, .id = "draw_col")
  # ecoregion scale
    ecoregion_ratios_list <- readRDS("./diversity_loss/ecoregion_ratios_list_100.rds")
    ecoregion_ratios_df <- bind_rows(lapply(ecoregion_ratios_list, function(draw) {
      bind_rows(draw, .id = "ecoregion")
    }), .id = "draw_col")
  # Pan-Colombia scale
    colombia_ratios_list <- readRDS("./diversity_loss/colombia_ratios_list_100.rds")
    colombia_ratios_df <- bind_rows(lapply(colombia_ratios_list, function(draw) {
      bind_rows(draw) %>% mutate(ecoregion = "Colombia")
    }), .id = "draw_col")

###########################
## Colombia-wide comparisonr
  # local scale 2km metrics
    mean(cell_ratios_df$med_logratio, na.rm = T)
    mean(exp(cell_ratios_df$med_logratio), na.rm = T)
    exp(mean(cell_ratios_df$med_logratio, na.rm = T))
    # ecoregions metrics
    mean(ecoregion_ratios_df$med_logratio, na.rm = T)
    mean(exp(ecoregion_ratios_df$med_logratio))
    exp(mean(ecoregion_ratios_df$med_logratio))
  # Pan-Colombia metrics  
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
  # IC pan_colombia
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
    
########
## maps 
  # log-ratios map
    ggplot(cell_ratios_df, aes(x = lon, y = lat, fill = avg_logratio)) +
      geom_raster() +
      scale_fill_viridis_c() +
      theme_bw() +
      labs(title = "Local scale diversity loss", fill = "Log ratio")
  # species richness map
    ggplot(cell_ratios_df, aes(x = lon, y = lat, fill = n)) +
      geom_raster() +
      scale_fill_viridis_c() +
      theme_bw() +
      labs(title = "Species richness", fill = "Number of species")

#########################################################
# ecoregions comparison across the posterior distribution 
  
  # Mean ecoregion ratios
    ecoregion_mean <- ecoregion_ratios_df %>%
      group_by(ecoregion) %>%
      summarise(med_logratio_mean = round(mean(med_logratio, na.rm = T),2),
                CI_5 = round(quantile(med_logratio, probs = 0.05, na.rm = TRUE),2),
                CI_95 = round(quantile(med_logratio, probs = 0.95, na.rm = TRUE),2))
    
  # 25th ecoregion ratios
    ecoregion_p25 <- ecoregion_ratios_df %>%
      group_by(ecoregion) %>%
      summarise(p25_logratio_mean = round(exp(mean(p_25_logratio, na.rm = T)),2),
                CI_5 = round(quantile(p_25_logratio, probs = 0.05, na.rm = TRUE),2),
                CI_95 = round(quantile(p_25_logratio, probs = 0.95, na.rm = TRUE),2))
  
  # 75th ecoregion ratios
    ecoregion_p75 <- ecoregion_ratios_df %>%
      group_by(ecoregion) %>%
      summarise(p_75_logratio_mean = round(mean(p_75_logratio, na.rm = T),2),
                CI_5 = round(quantile(p_75_logratio, probs = 0.05, na.rm = TRUE),2),
                CI_95 = round(quantile(p_75_logratio, probs = 0.95, na.rm = TRUE),2))
    
  # Merge ecoregion ratios with Pan-Colombia ratios
  ratios_df <- bind_rows(ecoregion_ratios_df, colombia_ratios_df)
  ratios_df$ecoregion <- gsub(" forests", "", ratios_df$ecoregion)
  ratios_df$ecoregion <- recode(ratios_df$ecoregion,
                              "Apure-Villavicencio dry"         = "Villavicencio dry",
                              "Cauca Valley montane"            = "Cauca montane",
                              "Cordillera Oriental montane"     = "EC montane",
                              "Eastern Cordillera Real montane" = "CC montane",
                              "Magdalena Valley dry"            = "Magdalena dry",
                              "Magdalena Valley montane"        = "Magdalena Valley",
                              "Northern Andean páramo"          = "Andean páramo",
                              "Northwest Andean montane"        = "WC montane")

  # boxplot ecoregion log ratios
    p1 <- ggplot(ratios_df, aes(x = reorder(ecoregion, p_25_logratio), y = p_25_logratio, fill = ecoregion)) +
      geom_boxplot(color = "black", alpha = 0.7) +
      scale_fill_viridis_d() +
      labs(title = "25th percentile", x = NULL, y = "Logratio") +
      theme_bw() +
      scale_y_continuous(limits=c(0,8)) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8))
    
    p2 <- ggplot(ratios_df, aes(x = reorder(ecoregion, avg_logratio), y = avg_logratio, fill = ecoregion)) +
      geom_boxplot(color = "black", alpha = 0.7) +
      scale_fill_viridis_d() +
      labs(title = "Mean", x = NULL, y = "Logratio") +
      theme_bw() +
      scale_y_continuous(limits=c(0,8)) +
      theme(legend.position = "none",axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8))
    
    p3 <- ggplot(ratios_df, aes(x = reorder(ecoregion, p_75_logratio), y = p_75_logratio, fill = ecoregion)) +
      geom_boxplot(color = "black", alpha = 0.7) +
      scale_fill_viridis_d() +
      labs(title = "75th percentile", x = NULL, y = "Logratio") +
      theme_bw() +
      scale_y_continuous(limits=c(0,8)) +
      theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8))
    
    grid.arrange(p1, p2, p3, ncol=3)

# Plot posterior distributions for metrics
    plot_25 <- ggplot(ratios_df, aes(x = p_25_logratio, y = reorder(ecoregion, -p_25_logratio), fill = ecoregion)) +
      geom_density_ridges(alpha = 0.5, scale = 1.5)  +
      scale_fill_viridis_d() +
      labs(title = "a", x = "sensitivity\n(N forest/N pasture)", y = "Ecoregion") +
      theme_bw() +
      theme(legend.position = "none",
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
      scale_x_continuous(limits=c(-1, 8)) +
      stat_density(geom = "line",position = "identity", aes(x = p_25_logratio), kernel = "gaussian")
    
    plot_mean <- ggplot(ratios_df, aes(x = med_logratio, y = reorder(ecoregion, -med_logratio), fill = ecoregion)) +
      geom_density_ridges(alpha = 0.5, scale = 1.5)  +
      scale_fill_viridis_d() +
      labs(title = "d", x = "sensitivity\n(N forest/N pasture)", y = "Ecoregion") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 9),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 10, face = "bold")) +
      scale_x_continuous(limits=c(-1, 8)) +
      stat_density(geom = "line", position = "identity", aes(x = med_logratio), kernel = "gaussian")
 
    plot_75 <- ggplot(ratios_df, aes(x = p_75_logratio, y = reorder(ecoregion, -p_75_logratio), fill = ecoregion)) +
      geom_density_ridges(alpha = 0.5, scale = 1.5)  +
      scale_fill_viridis_d() +
      labs(title = "g", x = "sensitivity (N forest/N pasture)", y = "Ecoregion") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 9),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 10, face = "bold")) +
      scale_x_continuous(limits=c(-1, 8)) +
      stat_density(geom = "line", position = "identity", aes(x = p_75_logratio), kernel = "gaussian")
    
    
    grid.arrange(plot_25, plot_mean, plot_75, ncol = 1)

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
      filter(ecoregion != "Colombia") %>% 
      arrange(med_logratio_dif_mean) %>%
      mutate(ecoregion_order = row_number())
    
    
    plot_mean_dif <- ggplot(relative_diff, aes(x = med_ratio_dif, y = reorder(ecoregion, med_ratio_dif), fill = ecoregion)) +
      geom_density_ridges(alpha = 0.5, scale = 1.5)  +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
      scale_fill_viridis_d() +
      labs(title = "e", x = "sensitivity\n(N forest/N pasture)", y = "Ecoregion") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 10),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 10, face = "bold")) +
      scale_x_continuous(limits = c(-1, 12), breaks = c(1, 4, 8, 12, 16, 20))+
      stat_density(geom = "line", position = "identity", aes(x = med_ratio_dif), kernel = "gaussian")
    
    
    plot_25_dif <- ggplot(relative_diff, aes(x = p_25_ratio_dif, y = reorder(ecoregion, p_25_ratio_dif), fill = ecoregion)) +
      geom_density_ridges(alpha = 0.5, scale = 1.5)  +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
      scale_fill_viridis_d() +
      labs(title = "b", x = "sensitivity\n(N forest/N pasture)", y = "Ecoregion") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 10),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 10, face = "bold")) +
      scale_x_continuous(limits = c(-1, 12), breaks = c(1, 4, 8, 12, 16, 20)) +
      stat_density(geom = "line",position = "identity", aes(x = p_25_ratio_dif), kernel = "gaussian")
    
    plot_75_dif <- ggplot(relative_diff, aes(x = p_75_ratio_dif, y = reorder(ecoregion, p_75_ratio_dif), fill = ecoregion)) +
      geom_density_ridges(alpha = 0.5, scale = 1.5)  +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
      scale_fill_viridis_d() +
      labs(title = "h", x = "Relative difference", y = "Ecoregion") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 10),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 10),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 10, face = "bold")) +
      scale_x_continuous(limits = c(-1, 12), breaks = c(1, 4, 8, 12, 16, 20)) +
      stat_density(geom = "line", position = "identity", aes(x = p_75_ratio_dif), kernel = "gaussian")
    
    grid.arrange(plot_25, plot_25_dif, plot_mean, 
                 plot_mean_dif, plot_75, plot_75_dif,
                 nrow = 3, widths =c(1.5,1))

#####################################################
## Plot Relative difference cummulative 
    
  # get summaries for metrics by ecoregions and Pan-Colombia sensitivities
    ratios_summary_df <- ratios_df %>%
      group_by(ecoregion) %>%
      summarise(
        avg_ratio = mean(avg_ratio, na.rm = TRUE),
        logratio_5 = quantile(med_logratio, probs = 0.05, na.rm = TRUE),
        logratio_95 = quantile(med_logratio, probs = 0.95, na.rm = TRUE),
        avg_logratio = mean(avg_logratio, na.rm = TRUE),
        med_logratio = mean(med_logratio, na.rm = TRUE),
        p_25_logratio = mean(p_25_logratio, na.rm = TRUE),
        p_75_logratio = mean(p_75_logratio, na.rm = TRUE),
        .groups = "drop") %>%
      filter(ecoregion != "Colombia") %>% 
      arrange(med_logratio) %>%
      mutate(ratio_exp = exp(med_logratio),
             CI_L = exp(logratio_5),
             CI_U = exp(logratio_95),
             ecoregion_order = row_number())  
    
  pan_colombia_avg <- colombia_ratios_df %>%
    group_by(draw_col) %>%  
    summarise(
      avg_logratio  = mean(avg_logratio, na.rm = TRUE),
      med_logratio  = mean(med_logratio, na.rm = TRUE),
      p_25_logratio  = mean(p_25_logratio, na.rm = TRUE),
      p_75_logratio  = mean(p_75_logratio, na.rm = TRUE),
      .groups = "drop")
  
  # select draws
  draws <- paste0("abun__draw_", 300 * c(1: 10))
  
  # random sequences
  set.seed(123)
  sequences <- replicate(1000, sample(unique(ecoregion_ratios_df$ecoregion)), simplify = FALSE)
  
  # loop for postrior distribution metrics
  
  results_all <- list()
  
  for (draw in draws) {
    cat("draw:", draw, "\n")
    
    col_draw <- pan_colombia_avg[pan_colombia_avg$draw_col == draw, ]
    ecoregion_data <- ecoregion_ratios_df[ecoregion_ratios_df$draw_col == draw, ]
    
    for (j in seq_along(sequences)) {
      cat("sequence", j, "\n")
      sequence <- sequences[[j]]
      
      metrics <- list(
        med_logratio = list(ref = col_draw$med_logratio, eco = ecoregion_data$med_logratio, metric = "median"),
        p_25_logratio = list(ref = col_draw$p_25_logratio, eco = ecoregion_data$p_25_logratio, metric = "p25"),
        p_75_logratio = list(ref = col_draw$p_75_logratio, eco = ecoregion_data$p_75_logratio, metric = "p75")
      )
      
      for (m in names(metrics)) {
        sens <- exp(metrics[[m]]$eco[match(sequence, ecoregion_data$ecoregion)])
        acc_sens <- sapply(seq_along(sens), function(i) mean(sens[1:i]))
        acc_log_sens <- log(acc_sens)
        
        rel_diff <- exp(metrics[[m]]$ref - acc_log_sens)
        
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
  
  # from list to dataframe
  results_df <- bind_rows(results_all)
  
  # summary by metrics
  summary_df <- results_df %>%
    group_by(ecoregion, metric) %>%
    summarise(
      mean_diff = mean(relative_difference),
      p_05 = quantile(relative_difference, probs = 0.05),
      p_95 = quantile(relative_difference, probs = 0.95),
      .groups = "drop")
  
  # Plotting metrics
  med_cum <- summary_df %>% filter(metric == "median")
  p25_cum <- summary_df %>% filter(metric == "p25")
  p75_cum <- summary_df %>% filter(metric == "p75")
  
  plot_p25_cum <- ggplot(p25_cum, aes(x = ecoregion, y = mean_diff)) +
    geom_point(size = 2.5, color = "black") +  # Puntos para la media
    geom_errorbar(aes(ymin = p_05, ymax = p_95), width = 0, color = "black", linewidth = 0.6) +  # Barras de error
    geom_ribbon(aes(ymin = 0.95, ymax = 1.05), fill = "gray", alpha = 0.5) +  # Banda gris
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Línea de referencia en 1
    scale_x_continuous(breaks = seq(1, 13, by = 2)) +  # Escala del eje X
    scale_y_continuous(limits=c(0.5, 6), breaks = seq(1, 6, by = 1)) +
    labs(title = "c", x = NULL, y = "Relative difference") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10))
  
  plot_med_cum <- ggplot(med_cum, aes(x = ecoregion, y = mean_diff)) +
    geom_point(size = 2.5, color = "black") +  # Puntos para la media
    geom_errorbar(aes(ymin = p_05, ymax = p_95), width = 0, color = "black", linewidth = 0.6) +  # Barras de error
    geom_ribbon(aes(ymin = 0.95, ymax = 1.05), fill = "gray", alpha = 0.5) +  # Banda gris
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Línea de referencia en 1
    scale_x_continuous(breaks = seq(1, 13, by = 2)) +  # Escala del eje X
    scale_y_continuous(limits=c(0.5, 8), breaks = seq(1, 8, by = 1)) +
    labs(title = "f", x = NULL, y = "Relative difference") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10))
  
  plot_p75_cum <- ggplot(p75_cum, aes(x = ecoregion, y = mean_diff)) +
    geom_point(size = 2.5, color = "black") +  # Puntos para la media
    geom_errorbar(aes(ymin = p_05, ymax = p_95),  width = 0, color = "black", linewidth = 0.6) +  # Barras de error
    geom_ribbon(aes(ymin = 0.95, ymax = 1.05), fill = "gray", alpha = 0.5) +  # Banda gris
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +  # Línea de referencia en 1
    scale_x_continuous(breaks = seq(1, 13, by = 2)) +  # Escala del eje X
    scale_y_continuous(limits=c(0.5, 9), breaks = seq(1, 8, by = 1)) +
    labs(title ="i", x = expression(N[regions]), y = "Relative difference") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size  = 14, face = "italic"),
      axis.title.y = element_text(size  = 10))
  
  fig_4 <- grid.arrange(plot_25, plot_25_dif, plot_p25_cum,
               plot_mean, plot_mean_dif, plot_med_cum,
               plot_75, plot_75_dif, plot_p75_cum,
               nrow = 3, widths =c(1.5,1, 1))
  
  ggsave("./fig_2.jpeg", plot = fig_2, width = 8.5, height = 9, units = "in",        
         dpi = 300, device = "jpeg")
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
 
    library(rnaturalearth)
    library(rnaturalearthdata)
    

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
    