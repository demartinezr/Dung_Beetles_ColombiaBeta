# ------------------------------------------------------------------------------
# Spatial biodiversity loss metrics across scales
#
# This script defines functions used to quantify biodiversity responses to
# forest–pasture conversion across spatial scales, following the procedures
# described in the Methods sections:
#
# -"Spatial scaling, abundance-based sensitivity metrics, and species inclusion criteria"
# -"Quantifying biodiversity loss across spatial scales using sensitivity to land-use conversion"
#
# The functions operate on predicted species abundances under forest and
# pasture scenarios and compute community-level sensitivity metrics at
# multiple spatial scales: local grid cells (2 × 2 km), ecoregions, and
# the near-national extent.
#
# Two complementary biodiversity indicators are implemented:
#
# 1) Abundance-based sensitivity
#    Median log-ratio of predicted abundances between forest and pasture
#    across species, representing the response of the median species in
#    a community.
#
# 2) Richness-based response
#    Proportion of species classified as declining ("losers"), defined as
#    species with a forest-to-pasture abundance ratio > 1.
#
# Species inclusion in each assemblage is controlled through minimum
# predicted abundance thresholds (cutoffs). The analyses support the
# thresholds used in the manuscript (1, 3, and 10 individuals), allowing
# sensitivity analyses of species pool definitions.
#
# The functions defined here are called by downstream scripts that compute
# biodiversity loss metrics from posterior predictions of the hierarchical
# abundance model.
# ------------------------------------------------------------------------------
setwd("C:/Users/PC/Dropbox/CO_DBdata")

library(dplyr)
library(data.table)
# ------------------------------------------------------------------------------
# avg_ratio
#
# Computes species-level forest-pasture abundance ratios and summarizes
# their distribution for a single spatial unit. Species predicted below the
# inclusion threshold are excluded from the assemblage. The function returns
# community-level summaries including the mean ratio, mean log-ratio,
# median log-ratio, and interquartile range of log-ratios across species.

avg_ratio <- function(fa_pa_i,  cutoff_use){
  n_col <- length(fa_pa_i)/2
  fa_i <- fa_pa_i[1:n_col]
  pa_i <- fa_pa_i[(n_col+1):(2*n_col)]
  incl <- ((fa_i) >= cutoff_use) | ((pa_i) >= cutoff_use)
  fa_inc <- fa_i[incl]
  pa_inc <- pa_i[incl]
  ratio <- fa_inc / pa_inc
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
# ------------------------------------------------------------------------------
# get_avg_cell_ratios
#
# Applies avg_ratio to each grid cell (2 × 2 km resolution) to calculate
# abundance-based sensitivity metrics at the local scale. For each cell,
# predicted forest and pasture abundances are compared across species and
# summarized as community-level statistics of the forest-pasture
# abundance ratios. Species inclusion is controlled by either absolute
# or relative abundance thresholds.

get_avg_cell_ratios <- function(forest, pasture, cutoff_type, cutoff){
  if(!all.equal(dim(forest), dim(pasture))){stop("forest_abun and pasture_abun have different dimensions")}
  if(cutoff <= 0){stop("cutoff must be greater than zero and less than one")}
  if(cutoff >= 12){stop("cutoff must be greater than zero and less than one")}
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
# ------------------------------------------------------------------------------
# get_regional_ratios
#
# Calculates abundance-based sensitivity metrics at broader spatial scales
# by aggregating predicted species abundances across multiple grid cells.
# When cell positions are provided, calculations are restricted to the
# specified spatial subset (e.g., an ecoregion). If no subset is provided,
# abundances are aggregated across all grid cells to obtain near-national
# estimates. Ratios are calculated from summed abundances and summarized
# using the same statistics applied at the local scale.

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
  log_ratio[is.infinite(log_ratio)] <- NA
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
# ------------------------------------------------------------------------------
# get_sample_percent_decline
#
# Computes richness-based biodiversity responses by identifying species
# that decline under pasture relative to forest. Species are classified as
# declining when the forest-to-pasture abundance ratio exceeds 1. The
# function returns the proportion of declining species at two levels:
# (1) pointwise estimates for individual grid cells, and
# (2) a pooled estimate based on aggregated abundances across the
# selected spatial scales.

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

