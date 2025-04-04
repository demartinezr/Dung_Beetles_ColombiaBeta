setwd("C:/Users/PC/Dropbox/CO_DBdata")

library(dplyr)
library(data.table)

# functions

# Function to compute the average ratio and log-ratio for a given set of presence-absence values.
avg_ratio <- function(fa_pa_i,  cutoff_use){
  n_col <- length(fa_pa_i)/2
  fa_i <- fa_pa_i[1:n_col]
  pa_i <- fa_pa_i[(n_col+1):(2*n_col)]
  incl <- ((fa_i) >= cutoff_use) | ((pa_i) >= cutoff_use)
  fa_inc <- fa_i[incl]
  pa_inc <- pa_i[incl]
  ratio <- fa_inc / pa_inc
  ratio[is.infinite(ratio)] <- 1
  ratio[is.nan(ratio)] <- NA
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
  if (cutoff <= 0) {stop("cutoff must be greater than zero and less than one")}
  if (cutoff >= 10) {stop("cutoff must be greater than zero and less than one")}
  
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