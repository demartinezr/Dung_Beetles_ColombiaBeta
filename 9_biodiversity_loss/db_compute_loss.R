setwd("C:/Users/PC/Dropbox/CO_DBdata")

library(dplyr)
library(data.table)

# functions

# Function to compute the average ratio and log-ratio for a given set of abundance values.
avg_ratio <- function(fa_pa_i,  cutoff_use){
  n_col <- length(fa_pa_i)/2
  fa_i <- fa_pa_i[1:n_col]
  pa_i <- fa_pa_i[(n_col+1):(2*n_col)]
  incl <- ((fa_i) >= cutoff_use) | ((pa_i) >= cutoff_use)
  fa_inc <- fa_i[incl]
  pa_inc <- pa_i[incl]
  ratio <- fa_inc / pa_inc
#  ratio[is.infinite(ratio)] <- 1
#  ratio[is.nan(ratio)] <- NA
  log_ratio <- log(ratio)
  log_ratio[is.infinite(log_ratio)] <- NA
  a_r <- mean(ratio, na.rm = TRUE)
  a_l <- mean(log_ratio, na.rm = TRUE)
  m_l <- median(log_ratio, na.rm = TRUE)
  p_25 <- quantile(log_ratio, probs = 0.25, na.rm = TRUE)
  p_75 <- quantile(log_ratio, probs = 0.75, na.rm = TRUE)
  n <- sum(incl)
  return(list(avg_ratio = a_r, avg_logratio = a_l, med_logratio = m_l, p_25_logratio = p_25, p_75_logratio = p_75, n = n))
}

get_avg_cell_ratios <- function(forest, pasture, cutoff_type, cutoff){
  if(!all.equal(dim(forest), dim(pasture))){stop("forest_abun and pasture_abun have different dimensions")}
  if(cutoff <= 0){stop("cutoff must be greater than zero and less than one")}
  if(cutoff >= 10){stop("cutoff must be greater than zero and less than one")}
  fa <- forest[, 4:ncol(forest)]
  pa <- pasture[, 4:ncol(pasture)]
  fa_pa <- cbind(fa, pa)
  
  if(cutoff_type == "relative"){
    ap <- rbind(fa, pa)
    maxabun <- apply(ap, 2, max)
    cutoff_use <- cutoff*maxabun
  }else if(cutoff_type == "absolute"){
    cutoff_use <- rep(cutoff, ncol(fa))
  }else{stop("cutoff type must be one of 'relative' or 'absolute'")}
  
  aln <- apply(fa_pa, 1, avg_ratio, cutoff_use = cutoff_use)
  output <- as.data.frame(do.call(rbind, aln))
  output$ecoregions <- as.factor(forest$ecoregions)
  output$avg_ratio <- as.numeric(output$avg_ratio)
  output$avg_logratio <- as.numeric(output$avg_logratio)
  output$med_logratio <- as.numeric(output$med_logratio)
  output$p_25_logratio <- as.numeric(output$p_25_logratio)
  output$p_75_logratio <- as.numeric(output$p_75_logratio)
  output$n <- as.integer(output$n)
  output$lon <- as.numeric(forest$lon)
  output$lat <- as.numeric(forest$lat)
    return(output)
}

get_regional_ratios <- function(forest, pasture,  cutoff_type, cutoff, cell_positions = NULL){
  if (!all.equal(dim(forest), dim(pasture))) {stop("forest_abun and pasture_abun have different dimensions")}
  if (cutoff <= 0) {stop("cutoff must be greater than zero and less than 0")}
  if (cutoff >= 12) {stop("cutoff must be greater than zero and less than 12")}
  
  fa <- forest[, 4:ncol(forest)]
  pa <- pasture[, 4:ncol(pasture)]
  
  if (cutoff_type == "relative") {
    ap <- rbind(fa, pa)
    maxabun <- apply(ap, 2, max)
    cutoff_use <- cutoff * maxabun
  } else if (cutoff_type == "absolute") {
    cutoff_use <- rep(cutoff, ncol(fa))
    }
 else {
      stop("cutoff type must be one of 'relative' or 'absolute'")
    }
  
  if (!is.null(cell_positions)) {
    fa <- fa[cell_positions, , drop = FALSE]
    pa <- pa[cell_positions, , drop = FALSE]
  }
  fa_max <- apply(fa, 2, max)
  pa_max <- apply(pa, 2, max)
  incl <- ((fa_max) >= cutoff_use) | ((pa_max) >= cutoff_use)
  
  fa_occ <- colSums(fa)[incl]
  pa_occ <- colSums(pa)[incl]
  ratio <- fa_occ / pa_occ
  ratio[is.infinite(ratio)] <- 1
  ratio[is.nan(ratio)] <- NA
  log_ratio <- log(ratio)
  log_ratio[is.infinite(log_ratio)] <- NA  # E -Inf
  avg_ratio <- mean(ratio, na.rm = TRUE)
  avg_logratio <- mean(log_ratio, na.rm = TRUE)
  med_logratio <- median(log_ratio, na.rm = TRUE)
  p_25_logratio <- quantile(log_ratio, probs = 0.25, na.rm = TRUE)
  p_75_logratio <- quantile(log_ratio, probs = 0.75, na.rm = TRUE)
  return(
    list(
      avg_ratio = avg_ratio,
      avg_logratio = avg_logratio,
      med_logratio = med_logratio,
      p_25_logratio = p_25_logratio,
      p_75_logratio = p_75_logratio,
      n = sum(incl)
    )
  )
}

get_sample_percent_decline <- function(forest, pasture, cutoff, cell_positions = NULL) {
  if(!all.equal(dim(forest), dim(pasture))){stop("forest and pasture have different dimensions")}
  if(cutoff < 1){stop("cutoff must be at least one")}
  fs <- as.data.frame(forest[, 4:ncol(forest), drop = FALSE])
  ps <-  as.data.frame(pasture[, 4:ncol(pasture), drop = FALSE])
  if(!is.null(cell_positions)){
    fs <- fs[cell_positions, , drop = FALSE]
    ps <- ps[cell_positions, , drop = FALSE]
  }
  f_plus_p <- fs + ps
  fs <- fs[, colSums(f_plus_p) >= cutoff, drop = FALSE]
  ps <- ps[, colSums(f_plus_p) >= cutoff, drop = FALSE]
  ratios_pointwise <- fs/ps
  decline_frac_pointwise <- apply(ratios_pointwise, 1, function(x){sum(x > 1, na.rm = T)/sum(!is.na(x))})
  ratios_total <- colSums(fs)/colSums(ps)
  nsp <- ncol(fs)
  decline_frac_total <- sum(ratios_total > 1)/nsp
  return(list(pointwise = decline_frac_pointwise, total = decline_frac_total, nsp = nsp))
}

################################################################################
# hexagons approach

# Generate wide-format species abundance data for forest and pasture across multiple draws
generate_forest_pasture_lists <- function(ecoregions_all, draws) {
  forest_list <- vector("list", length(draws))
  pasture_list <- vector("list", length(draws))
  names(forest_list) <- names(pasture_list) <- paste0("draw_", draws)
  
  for (draw in draws) {
    cat("draw", draw, "\n")
    draw_col <- paste0("abun__draw_", draw)
    
    ecoregions_tmp <- ecoregions_all %>%
      select(ecoregions, lon, lat, scientificName, pasture, all_of(draw_col)) %>%
      mutate(cell_id = row_number()) %>%
      select(cell_id, everything())
    
    row.names(ecoregions_tmp) <- NULL
    
    forest_tmp <- ecoregions_tmp %>% filter(pasture == 0) %>% select(-pasture)
    
    pasture_tmp <- ecoregions_tmp %>% filter(pasture == 1) %>% select(-pasture)
    
    forest_list[[paste0("forest_", draw, "_wide")]] <- dcast(as.data.table(forest_tmp),
                                                             ecoregions + lon + lat ~ scientificName,
                                                             value.var = draw_col, fill = 0)
    
    pasture_list[[paste0("pasture_", draw, "_wide")]] <- dcast(as.data.table(pasture_tmp),
                                                               ecoregions + lon + lat ~ scientificName,
                                                               value.var = draw_col, fill = 0)
  }
  
  return(list(forest_list = forest_list, pasture_list = pasture_list))
}

# Function to compute diversity metrics at the cell level
process_cells <- function(forest_data, pasture_data, cell_ratios, unique_cells, cell_ids) {
  n_cells <- length(unique_cells)
  n_points <- length(cell_ids)
  
  # Initialize result vectors
  regional_logratios <- numeric(n_cells)
  regional_richness <- numeric(n_cells)
  pointwise_logratios <- numeric(n_cells)
  pointwise_richness <- numeric(n_cells)
  beta_diversity <- numeric(n_cells)
  cell_counts <- numeric(n_cells)
  
  # Expanded vectors (same length as number of points)
  expanded_logratios <- rep(NA, n_points)
  expanded_pointwise_logratios <- rep(NA, n_points)
  expanded_beta <- rep(NA, n_points)
  
  for (i in seq_along(unique_cells)) {
    cell_value <- unique_cells[i]
    cell_positions <- which(cell_ids == cell_value)
    
    regional <- get_regional_ratios(forest_data, pasture_data, cutoff_type = "absolute", cutoff = 1, cell_positions = cell_positions)
    
    regional_logratios[i] <- regional$med_logratio
    regional_richness[i] <- regional$n
    expanded_logratios[cell_positions] <- regional_logratios[i]
    
    # Filtro correcto solo en las posiciones de esta celda
    valid_positions <- cell_positions[!is.na(cell_ratios$avg_logratio[cell_positions])]
    pointwise_logratios[i] <- mean(cell_ratios$med_logratio[cell_positions], na.rm = TRUE)
    pointwise_richness[i] <- mean(cell_ratios$n[valid_positions])
    beta_diversity[i] <- regional_richness[i] / pointwise_richness[i]
    
    expanded_pointwise_logratios[cell_positions] <- pointwise_logratios[i]
    expanded_beta[cell_positions] <- beta_diversity[i]
    cell_counts[i] <- length(cell_positions)
  }
  
  return(list(
    expanded_logratios = expanded_logratios,
    expanded_pointwise_logratios = expanded_pointwise_logratios,
    expanded_beta = expanded_beta
  ))
}

# Function to create raster layers from diversity metrics
create_rasters <- function(coords, logratios, pointwise, beta) {
  raster_regional <- raster::rasterFromXYZ(cbind(coords, logratios))
  raster_pointwise <- raster::rasterFromXYZ(cbind(coords, pointwise))
  raster_beta <- raster::rasterFromXYZ(cbind(coords, beta))
  raster_diff <- raster_regional - raster_pointwise
  
  return(list(
    raster_regional = raster_regional,
    raster_pointwise = raster_pointwise,
    raster_beta = raster_beta,
    raster_difference = raster_diff
  ))
}

# Function to handle one draw-resolution combination
process_draw_resolution <- function(draw, res_index, resolutions, dggs_list, cell_list, unique_cells_list,
                                    forest_list, pasture_list, coords) {
  resolution <- resolutions[res_index]
  draw_id <- paste0("draw_", draw)
  
  forest_data <- forest_list[[paste0("forest_", draw, "_wide")]]
  pasture_data <- pasture_list[[paste0("pasture_", draw, "_wide")]]
  
  cell_ratios <- get_avg_cell_ratios(forest_data, pasture_data, cutoff_type = "absolute", cutoff = 1)
  cell_ids <- cell_list[[res_index]]
  unique_cells <- unique_cells_list[[res_index]]
  
  cell_results <- process_cells(forest_data, pasture_data, cell_ratios, unique_cells, cell_ids)
  
  rasters <- create_rasters(
    coords = forest_data[, c("lon", "lat")],
    logratios = cell_results$expanded_logratios,
    pointwise = cell_results$expanded_pointwise_logratios,
    beta = cell_results$expanded_beta
  )
  
  return(rasters)
}

# convert raster to data frame con renombramiento manual
raster_to_df <- function(raster, resolution, draw, type, colname) {
  df <- as.data.frame(raster, xy = TRUE)
  
  df %>%
    rename(logratio = all_of(colname)) %>%
    mutate(res = paste0(resolution),
           draw = draw,
           type = type)
}