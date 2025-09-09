setwd("C:/Users/PC/Dropbox/CO_DBdata")
ecoregions_predictions <- readRDS("./Analysis/mean_abundance/ecoregions_predictions.rds")

library(raster)
library(purrr)
library(dplyr)
library(data.table)

AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

# ecoregions names in each sf dataframe
ecoregions_predictions <- map2(
  ecoregions_predictions, 
  names(ecoregions_predictions), 
  ~ mutate(.x, ecoregions = .y)
)
# from list to sf dataframe
ecoregions_all <- reduce(ecoregions_predictions, rbind)
# filter draws
ecoregions_all <- as.data.frame(readRDS("./ecoregions_all.rds"))

ecoregions_900 <- ecoregions_all %>%
  select(ecoregions, lon, lat, scientificName, pasture, abun__draw_1800) %>%
  mutate(cell_id = row_number())%>%
  select(cell_id, everything()) 
  row.names(ecoregions_900) <- NULL
 
  forest_900 <- ecoregions_900 %>%
    filter(pasture == 0)  %>%
    select(-pasture)
  
  forest_900_dt <- as.data.table(forest_900)
  forest_900_wide <- dcast(forest_900_dt, ecoregions + lon + lat ~ scientificName, 
                            value.var = "abun__draw_1800", fill = 0)
  
  pasture_900 <- ecoregions_900 %>%
    filter(pasture == 1)  %>%
    select(-pasture)
  
  pasture_900_dt <- as.data.table(pasture_900)
  pasture_900_wide <- dcast(pasture_900_dt, ecoregions + lon + lat ~ scientificName, 
                              value.var = "abun__draw_1800", fill = 0)
  
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
    fa <- forest[, 3:ncol(forest)]
    pa <- pasture[, 3:ncol(pasture)]
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
  
cell_ratios <- get_avg_cell_ratios(forest_900_wide, pasture_900_wide, cutoff_type="absolute", cutoff=1)
  
  get_regional_ratios <- function(forest, pasture, cutoff_type, cutoff, cell_positions = NULL){
    if(!all.equal(dim(forest), dim(pasture))){stop("forest_probs and pasture_probs have different dimensions")}
    if(cutoff <= 0){stop("cutoff must be greater than zero and less than one")}
    if(cutoff >= 10){stop("cutoff must be greater than zero and less than one")}
    fa <- forest[, 4:ncol(forest)]
    pa <- pasture[, 4:ncol(pasture)]
   
    if(cutoff_type == "relative"){
      ap <- rbind(fa, pa)
      maxprobs <- apply(ap, 2, max)
      cutoff_use <- cutoff*maxprobs
    }else if(cutoff_type == "absolute"){
      cutoff_use <- rep(cutoff, ncol(fa))
    }else{stop("cutoff type must be one of 'relative' or 'absolute'")}
   
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
    log_ratio[is.infinite(log_ratio)] <- NA  # Evita -Inf
    avg_ratio <- mean(ratio, na.rm = TRUE)
    avg_logratio <- mean(log_ratio, na.rm = TRUE)
    med_logratio <- median(log_ratio, na.rm = TRUE)
    return(list(avg_ratio = avg_ratio, avg_logratio=avg_logratio, med_logratio=med_logratio, n = sum(incl)))
  }
  
colombia_ratio <- get_regional_ratios(forest_900_wide, pasture_900_wide, cutoff_type="absolute", cutoff=1, cell_positions = NULL)

##### Colombia-wide comparison #####
mean(cell_ratios$med_logratio, na.rm = T)
mean(exp(cell_ratios$med_logratio), na.rm = T)
exp(mean(cell_ratios$med_logratio, na.rm = T))
colombia_ratio
exp(colombia_ratio$med_logratio)


# Map the local averages
cell_logratios_df <- data.frame(x=forest_900_wide$lon, y = forest_900_wide$lat, logratio=cell_ratios$avg_logratio)
cell_logratios_raster <- rasterFromXYZ(cell_logratios_df)
par(mar=c(2,2,1,1))
plot(cell_logratios_raster)

# Map species richness
cell_richness_df <- data.frame(x=forest_900_wide$lon, y=forest_900_wide$lat, richness = cell_ratios$n)
cell_richness_raster <- raster::rasterFromXYZ(cell_richness_df)
raster::plot(cell_richness_raster)
