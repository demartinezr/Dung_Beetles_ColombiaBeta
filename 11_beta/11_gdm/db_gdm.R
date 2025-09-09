# Script to fit gdms to raw data and posterior matrix
library(dplyr)
library(reticulate)
library(gdm)
library(sf)
library(gtools)
library(raster)
library(betapart)
library(rgee)

setwd("C:/Users/PC/Dropbox/CO_DBdata")
##### Formating covariates #####
# dataset with abundance by speices/point/day -> 923315 obs, 958 points, 243 species
db5 <- readRDS("abundance/db5_distance.RDS")
db5$elev_standard_squared <- db5$elev_standard^2
db5$subregion_species <- paste0(db5$subregion, "__", db5$scientificName)
db5$cluster_species <- paste0(db5$cluster, "__", db5$scientificName)
db5$distance_from_range_scaled2 <- as.vector(db5$distance_from_range_scaled2)
# sum abundance day by species/point. Remove day covariate -> 241948 obs, 958 points, 243 species
db6 <- db5 %>%  
  group_by(across(-c(day, p_d))) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")  
# remove palm and young forest 219452 obs, 870 points, 243 species
db6 <- db6[db6$habitat !="Sy",]
db6 <- db6[db6$habitat !="PALM",]

# Remove the improbable combinations for point and species range, 865 points, 243 species
db7 <- db6 |>
  filter(combinations ==1)

class(db7$pasture)
db7$pasture <- as.factor(db7$pasture)

db7_pts <- db7[!duplicated(db7$point), c("point", "lat", "lon", "elev_ALOS", "pasture")]

## use gee to get precipitation data, and save result
reticulate::use_condaenv("r-reticulate", required = TRUE)
## ee_Authenticate()
ee_Initialize()

 BioClim <- ee$Image('WORLDCLIM/V1/BIO')
 BioClim_precip <- BioClim$select('bio12')
 geompts <- sapply(1:nrow(db7_pts),function(x)ee$Geometry$Point(c(db7_pts$lon[x],db7_pts$lat[x])))
 geompts <- ee$FeatureCollection(c(unlist(geompts)))
 pts_precip <- BioClim$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
 BioClim_precip <- sapply(c(1:length(pts_precip$features)),function(x)pts_precip$features[[x]]$properties$bio12)
 db7_pts$precip <- BioClim_precip
 CHIRP <- raster("F:/Capas/America/weather/SA_MAP_1km_CHIRP/SA_MAP_1km_CHIRP.tif")
 db7_pts$precip_ceccherini <- extract(CHIRP, cbind(db7_pts$lon, db7_pts$lat))
 IDEAM <- raster("F:/Capas/America/weather/Escenario_Precipitacion_1976_2005/ECC_Prcp_GeoTiff_2011_2040/ECC_Prcp_1976_2005_100K_2015.tif")
 db7_pts$precip_IDEAM <- extract(IDEAM, cbind(db7_pts$lon, db7_pts$lat))
 db7_pts[db7_pts$precip_ceccherini<1000,]
 db7_pts[db7_pts$precip<1000,]
 db7_pts[grep("CC", db7_pts$point),]
 saveRDS(db7_pts, "./gdm/db7_pts_w_precip.RDS")

db7_pts <- readRDS("./gdm/db7_pts_w_precip.RDS")

## get biogeographic regions for each point
#source("F:/repositorio/Dung_Beetles_ColombiaBeta/3_Geographic_range/Topographic_units/unidades_topograficas.R")

#bp2 <- st_as_sf(db7_pts, coords = c("lon", "lat"))
#st_crs(bp2) <- 4326
#east <- st_union(catatumbo, amazon_orinoco)
#eAndes <- st_union(east, magdalena_east)
#eAndes <- st_make_valid(eAndes)
#cAndes <- st_union(magdalena_west, cauca_east)
#wAndes <- st_union(pacific, cauca_west)
#sf_use_s2(FALSE)
#all_regions <- st_union(eAndes, st_union(cAndes, st_union(wAndes, st_union(snsm, pasto))))
#if(sum(!st_covered_by(bp2, all_regions, sparse=F)) == 0){print("All points assigned to a region :)")}
#db7_pts$east <- st_covered_by(bp2, east, sparse = F)
#db7_pts$range <- "0"
#db7_pts$range[st_covered_by(bp2, eAndes, sparse = F) | st_covered_by(bp2, pasto, sparse = F)] <- "east_south"
#db7_pts$range[st_covered_by(bp2, cAndes, sparse = F)] <- "central"
#db7_pts$range[st_covered_by(bp2, wAndes, sparse = F)] <- "west"
#db7_pts$range[st_covered_by(bp2, snsm, sparse = F)] <- "snsm"
#saveRDS(db7_pts, "./gdm/db7_pts.RDS")
db7_pts <- readRDS("./gdm/db7_pts.RDS")

db7 <- merge(
  db7,
  db7_pts[, c("point", "precip", "precip_ceccherini", "precip_IDEAM", "east", "range")],
  by = "point",
  all.x = TRUE
)

# load dataset and restrict to points with a detection
db8 <- db7[db7$abundance >= 1, ]

# Extract points and some covariates
db8_pts <- db8[!duplicated(db8$point), c("point", "lat", "lon", "elev_ALOS", "pasture")]

#BioClim <- ee$Image('WORLDCLIM/V1/BIO')
#BioClim_precip <- BioClim$select('bio12')
#geompts <- sapply(1:nrow(db8_pts),function(x)ee$Geometry$Point(c(db8_pts$lon[x],db8_pts$lat[x])))
#geompts <- ee$FeatureCollection(c(unlist(geompts)))
#pts_precip <- BioClim$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
#BioClim_precip <- sapply(c(1:length(pts_precip$features)),function(x)pts_precip$features[[x]]$properties$bio12)
#db8_pts$precip <- BioClim_precip
#CHIRP <- raster("F:/Capas/America/weather/SA_MAP_1km_CHIRP/SA_MAP_1km_CHIRP.tif")
#db8_pts$precip_ceccherini <- extract(CHIRP, cbind(db8_pts$lon, db8_pts$lat))
#IDEAM <- raster("F:/Capas/America/weather/Escenario_Precipitacion_1976_2005/ECC_Prcp_GeoTiff_2011_2040/ECC_Prcp_1976_2005_100K_2015.tif")
#db8_pts$precip_IDEAM <- extract(IDEAM, cbind(db8_pts$lon, db8_pts$lat))
#db8_pts[db8_pts$precip_ceccherini<1000,]
#db8_pts[db8_pts$precip<1000,]
#db8_pts[grep("CC", db8_pts$point),]
#saveRDS(db8_pts, "./gdm/db8_pts_w_precip.RDS")

 db8_pts <- readRDS("./gdm/db8_pts_w_precip.RDS")

# get biogeographic regions for each point
source("F:/repositorio/Dung_Beetles_ColombiaBeta/3_Geographic_range/Topographic_units/unidades_topograficas.R")

bp2 <- st_as_sf(db8_pts, coords = c("lon", "lat"))
st_crs(bp2) <- 4326
east <- st_union(catatumbo, amazon_orinoco)
eAndes <- st_union(east, magdalena_east)
eAndes <- st_make_valid(eAndes)
cAndes <- st_union(magdalena_west, cauca_east)
wAndes <- st_union(pacific, cauca_west)
sf_use_s2(FALSE)
all_regions <- st_union(eAndes, st_union(cAndes, st_union(wAndes, st_union(snsm, pasto))))
if(sum(!st_covered_by(bp2, all_regions, sparse=F)) == 0){print("All points assigned to a region :)")}
db8_pts$east <- st_covered_by(bp2, east, sparse = F)
db8_pts$range <- "0"
db8_pts$range[st_covered_by(bp2, eAndes, sparse = F) | st_covered_by(bp2, pasto, sparse = F)] <- "east_south"
db8_pts$range[st_covered_by(bp2, cAndes, sparse = F)] <- "central"
db8_pts$range[st_covered_by(bp2, wAndes, sparse = F)] <- "west"
db8_pts$range[st_covered_by(bp2, snsm, sparse = F)] <- "snsm"
#saveRDS(db8_pts, "./gdm/db8_pts.RDS")

# Separate to forest points and pasture points
db8_pts <- readRDS("./gdm/db8_pts.RDS")
db8_pts_f <- db8_pts_f_stable <- db8_pts[db8_pts$pasture == 0, names(db8_pts) != "pasture"] # the "stable" version will be a version to come back to when we start modifying birdpts_f later on
db8_pts_p <- db8_pts_p_stable <- db8_pts[db8_pts$pasture == 1, names(db8_pts) != "pasture"]
db8_f <- db8[db8$pasture == 0, c("scientificName", "point", "lat", "lon")]
db8_p <- db8[db8$pasture == 1, c("scientificName", "point", "lat", "lon")]
db8_pts_p <- db8_pts_p[db8_pts_p$point %in% db8_p$point,] # exclude points with no detections in the raw pasture data (there are no such points in the forest data)
points_f <- db8_pts_f$point
points_p <- db8_pts_p$point

# For forest points:
#   distance matrix for montane barriers
montane_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
for(i in 1:length(points_f)){
  for(j in 1:length(points_f)){
    if(db8_pts_f$east[i] != db8_pts_f$east[j]){
      montane_matrix_f[i,j] <- (4100 - max(db8_pts_f$elev_ALOS[i], db8_pts_f$elev_ALOS[j]))
    }
  }
}
montane_barrier_f <- cbind(points_f, data.frame(montane_matrix_f))
names(montane_barrier_f)[1] <- "point"
#   distance matrix for valley barriers
valley_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
for(i in 1:length(points_f)){
  for(j in 1:length(points_f)){
    if(db8_pts_f$range[i] != db8_pts_f$range[j]){
      valley_matrix_f[i,j] <- (min(db8_pts_f$elev_ALOS[i], db8_pts_f$elev_ALOS[j]))
    }
  }
}
valley_barrier_f <- cbind(points_f, data.frame(valley_matrix_f))
names(valley_barrier_f)[1] <- "point"

# Pasture points
#   distance matrix for montane barriers
montane_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
for(i in 1:length(points_p)){
  for(j in 1:length(points_p)){
    if(db8_pts_p$east[i] != db8_pts_p$east[j]){
      montane_matrix_p[i,j] <- (4100 - max(db8_pts_p$elev_ALOS[i], db8_pts_p$elev_ALOS[j]))
    }
  }
}
montane_barrier_p <- cbind(points_p, data.frame(montane_matrix_p))
names(montane_barrier_p)[1] <- "point"
#   distance matrix for valley barriers
valley_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
for(i in 1:length(points_p)){
  for(j in 1:length(points_p)){
    if(db8_pts_p$range[i] != db8_pts_p$range[j]){
      valley_matrix_p[i,j] <- (min(db8_pts_p$elev_ALOS[i], db8_pts_p$elev_ALOS[j]))
    }
  }
}
valley_barrier_p <- cbind(points_p, data.frame(valley_matrix_p))
names(valley_barrier_p)[1] <- "point"


########################## Raw data gdm: sorensen ##############################

# Format site-pair tables for forest and pasture
gdmTab_f <- formatsitepair(db8_f, bioFormat=2, XColumn="lon", YColumn="lat",
                           sppColumn="scientificName", siteColumn="point", 
                           predData = as.data.frame(db8_pts_f[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")]),
                           distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
gdmTab_p <- formatsitepair(db8_p, bioFormat=2, XColumn="lon", YColumn="lat",
                           sppColumn="scientificName", siteColumn="point", 
                           predData = as.data.frame(db8_pts_p[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")]),
                           distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))



# Fit GDM for forest points
#   Observed data
forest_gdm_obs <- gdm(gdmTab_f, geo = T)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_f <- gdmTab_f
forest_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
  pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
  withweights_f$weights <- pair_weights_f
  forest_gdm_bb[[i]] <- gdm(withweights_f, geo = T)
}

# Fit GDM for pasture points
#   Observed data
pasture_gdm_obs <- gdm(gdmTab_p, geo = T)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_p <- gdmTab_p
pasture_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
  pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
  withweights_p$weights <- pair_weights_p
  pasture_gdm_bb[[i]] <- gdm(withweights_p, geo = T)
}

gdms_raw <- list(forest_obs = forest_gdm_obs, pasture_obs = pasture_gdm_obs, forest_bb = forest_gdm_bb, pasture_bb = pasture_gdm_bb)
saveRDS(gdms_raw, file = "C:/Users/PC/Dropbox/CO_DBdata/gdm/gdms_raw.RDS")

############################ Modeled data GDM #################################

##### do prediction #####
library(posterior)
library(brms)

db_traits <- readRDS("./abundance/db5_traits.RDS") |>
  filter(!duplicated(scientificName)) |>
  dplyr::select(scientificName, sp_elev_lower2, sp_elev_upper2, sp_region_amazon,
         sp_region_llanos, sp_region_caribbean, sp_region_snsm, sp_region_andes,
         sp_region_eastern, sp_slope_ECe, sp_slope_ECw, sp_slope_CCe, sp_slope_CCw,
         sp_slope_WCe, sp_slope_WCw, sp_slope_SNSM, nest_guild, diet_range, activity,
         bodysize, legratio)

db7_pts <- readRDS("./gdm/db7_pts.RDS")

db7 <- merge(
  db7,
  db7_pts[, c("point", "precip", "precip_ceccherini", "precip_IDEAM", "east", "range")],
  by = "point",
  all.x = TRUE
)

mod <- readRDS("./db_mod_abundance.rds")
mod_draws <- as_draws_df(mod)

species_preds <- cluster_sds <- phis <- list()
draw_ids <- 30 * c(1:100)
for(i in 1:nrow(db_traits)){
  print(i)
  species_preds[[i]] <- posterior_linpred(
    mod, 
    draw_ids = draw_ids,
    re_formula = ~ (1 + pasture + elev_standard + elev_standard_squared | scientificName) +
      (1 | subregion_species),
    newdata = db7 |>
      filter(scientificName == db_traits$scientificName[i]),
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian"
  )
  cluster_sds[[i]] <- mod_draws$sd_cluster_species__Intercept[draw_ids]
  phis[[i]] <- mod_draws$shape[draw_ids]
}

simulate_abun <- Vectorize(function(linpred, cluster_sd, phi) {
  out <- 0
  for(i in 1:16){
    lp <- linpred + rnorm(1, 0, cluster_sd)
    abun <- rnbinom(5, size = phi, mu = exp(lp))
    out <- out + mean(abun)/16
  }
  out
})

species_predictions <- list()
for(i in 1:nrow(db_traits)){
  print(i)
  pd_species <- db7 |>
    filter(scientificName == db_traits$scientificName[i])
  predictions <- as.data.frame(matrix(nrow = nrow(pd_species), ncol = length(draw_ids)))
  names(predictions) <- paste0("abun__draw_", draw_ids)
  for(j in 1:length(draw_ids)) {
    print(j)
    predictions[,j] <- 
      simulate_abun(species_preds[[i]][j, ], cluster_sds[[i]][j], phis[[i]][j])
  }
  species_predictions[[i]] <- cbind(pd_species, predictions)
}

names(species_predictions) <- db_traits$scientificName
# saveRDS(species_predictions, "./gdm/species_predictions_gdm.RDS")


species_predictions <- readRDS("./gdm/species_predictions_gdm.RDS")

a_rep_df <- bind_rows(species_predictions)

# Nombres de las columnas de abundancia
abun_cols <- grep("^abun__draw_", colnames(a_rep_df), value = TRUE)

# Crear lista con un dataframe por cada columna de abundancia
a_rep <- lapply(abun_cols, function(col_name) {
  df <- a_rep_df[, setdiff(colnames(a_rep_df), abun_cols)] # todas las columnas excepto abun__draw_XX
  df$a <- a_rep_df[[col_name]]  # nueva columna 'a' con la abundancia específica
  return(df)
})

# Asignar nombres a la lista
names(a_rep) <- abun_cols

forest_gdm_rep <- pasture_gdm_rep <- list()
for(k in 1:100){ # k indexes the posterior iteration
  print(k)
  db9 <- db7[a_rep[[k]]$a >= 1, ]
  db9_pts <- db9[!duplicated(db9$point), c("point", "lat", "lon", "elev_ALOS", "pasture", "precip","precip_ceccherini", "precip_IDEAM", "east", "range")]
  db9_f <- db9[db9$pasture == 0, c("scientificName", "point", "lat", "lon")]
  db9_pts_f <- db9_pts_f_stable <- db9_pts[db9_pts$pasture == 0, names(db9_pts) != "pasture"]
  db9_p <- db9[db9$pasture == 1, c("scientificName", "point", "lat", "lon")]
  db9_pts_p <- db9_pts_p_stable <- db9_pts[db9_pts$pasture == 1, names(db9_pts) != "pasture"]
  db9pts_f <- db9_pts_f_stable[db9_pts_f_stable$point %in% db9_f$point,] # exclude points with no detections
  db9pts_p <- db9_pts_p_stable[db9_pts_p_stable$point %in% db9_p$point,] # exclude points with no detections
  # forest
  points_f <- db9pts_f$point
  montane_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
  for(i in 1:length(points_f)){
    for(j in 1:length(points_f)){
      if(db9pts_f$east[i] != db9pts_f$east[j]){
        montane_matrix_f[i,j] <- (4100 - max(db9pts_f$elev_ALOS[i], db9pts_f$elev_ALOS[j]))
      }
    }
  }
  montane_barrier_f <- cbind(points_f, data.frame(montane_matrix_f))
  names(montane_barrier_f)[1] <- "point"

  valley_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
  for(i in 1:length(points_f)){
    for(j in 1:length(points_f)){
      if(db9pts_f$range[i] != db9pts_f$range[j]){
        valley_matrix_f[i,j] <- (min(db9pts_f$elev_ALOS[i], db9pts_f$elev_ALOS[j]))
      }
    }
  }
  valley_barrier_f <- cbind(points_f, data.frame(valley_matrix_f))
  names(valley_barrier_f)[1] <- "point"

  # pasture
  points_p <- db9pts_p$point
  montane_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
  for(i in 1:length(points_p)){
    for(j in 1:length(points_p)){
      if(db9pts_p$east[i] != db9pts_p$east[j]){
        montane_matrix_p[i,j] <- (4100 - max(db9pts_p$elev_ALOS[i], db9pts_p$elev_ALOS[j]))
      }
    }
  }
  montane_barrier_p <- cbind(points_p, data.frame(montane_matrix_p))
  names(montane_barrier_p)[1] <- "point"


  valley_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
  for(i in 1:length(points_p)){
    for(j in 1:length(points_p)){
      if(db9pts_p$range[i] != db9pts_p$range[j]){
        valley_matrix_p[i,j] <- (min(db9pts_p$elev_ALOS[i], db9pts_p$elev_ALOS[j]))
      }
    }
  }
  valley_barrier_p <- cbind(points_p, data.frame(valley_matrix_p))
  names(valley_barrier_p)[1] <- "point"

  gdmTab_f <- formatsitepair(db9_f, bioFormat=2, XColumn="lon", YColumn="lat",
                             sppColumn="scientificName", siteColumn="point", 
                             predData = as.data.frame(db9pts_f[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")]),
                             distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
  gdmTab_p <- formatsitepair(db9_p, bioFormat=2, XColumn="lon", YColumn="lat",
                             sppColumn="scientificName", siteColumn="point", 
                             predData = as.data.frame(db9pts_p[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")]),
                             distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))
  
  # Bootstrap replicate (one per posterior iteration)
  n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
  n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
  s_all <- unique(c(n1,n2))
  s1 <- match(n1, s_all)
  s2 <- match(n2, s_all)
  withweights_f <- gdmTab_f
  site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
  pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
  withweights_f$weights <- pair_weights_f
  
  n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
  n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
  s_all <- unique(c(n1,n2))
  s1 <- match(n1, s_all)
  s2 <- match(n2, s_all)
  withweights_p <- gdmTab_p
  site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
  pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
  withweights_p$weights <- pair_weights_p
  
  forest_gdm_rep[[k]] <- gdm(withweights_f, geo = T)
  pasture_gdm_rep[[k]] <- gdm(withweights_p, geo = T)
}

gdms_modeled <- list(forest_gdm_rep_bb = forest_gdm_rep, pasture_gdm_rep_bb = pasture_gdm_rep)
saveRDS(gdms_modeled, file = "./gdm/gdms_modeled_v6.RDS")

################################################################################
#################### Raw data gdm: turnover (Simpson) #########################
db8_pts <- readRDS("./gdm/db8_pts.RDS")
db8_pts_f <- db8_pts_f_stable <- db8_pts[db8_pts$pasture == 0, names(db8_pts) != "pasture"] # the "stable" version will be a version to come back to when we start modifying birdpts_f later on
db8_pts_p <- db8_pts_p_stable <- db8_pts[db8_pts$pasture == 1, names(db8_pts) != "pasture"]
db8_f <- db8[db8$pasture == 0, c("scientificName", "point", "lat", "lon")]
db8_p <- db8[db8$pasture == 1, c("scientificName", "point", "lat", "lon")]
db8_pts_p <- db8_pts_p[db8_pts_p$point %in% db8_p$point,] # exclude points with no detections in the raw pasture data (there are no such points in the forest data)
points_f <- db8_pts_f$point
points_p <- db8_pts_p$point

# Format site-pair tables for forest and pasture
db_spp_f <- unique(db8_f$scientificName)
db_matrix_f <- matrix(data=0, nrow=nrow(db8_pts_f), ncol = length(db_spp_f))
db_matrix_f[cbind(match(db8_f$point, db8_pts_f$point), match(db8_f$scientificName, db_spp_f))] <- 1
forest_betapart <- beta.pair(db_matrix_f)
db_format3_f <- cbind(db8_pts_f$point, as.matrix(forest_betapart$beta.sim))
colnames(db_format3_f)[1] <- "point"
db_format3_f[, 1] <- 1:nrow(db_format3_f)
db_format3_f <- as.data.frame(db_format3_f)
db_format3_f$point <- as.numeric(db_format3_f$point)

db_spp_p <- unique(db8_p$scientificName)
db_matrix_p <- matrix(data=0, nrow=nrow(db8_pts_p), ncol = length(db_spp_p))
db_matrix_p[cbind(match(db8_p$point, db8_pts_p$point), match(db8_p$scientificName, db_spp_p))] <- 1
pasture_betapart <- beta.pair(db_matrix_p)
db_format3_p <- cbind(db8_pts_p$point, as.matrix(pasture_betapart$beta.sim))
colnames(db_format3_p)[1] <- "point"
db_format3_p[, "point"] <- 1:nrow(db_format3_p)
db_format3_p <- as.data.frame(db_format3_p)
db_format3_p$point <- as.numeric(db_format3_p$point)

# For forest points:
#   distance matrix for montane barriers
montane_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
for(i in 1:length(points_f)){
  for(j in 1:length(points_f)){
    if(db8_pts_f$east[i] != db8_pts_f$east[j]){
      montane_matrix_f[i,j] <- (4100 - max(db8_pts_f$elev_ALOS[i], db8_pts_f$elev_ALOS[j]))
    }
  }
}
montane_barrier_f <- cbind(points_f, data.frame(montane_matrix_f))
names(montane_barrier_f)[1] <- "point"

#   distance matrix for valley barriers
valley_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
for(i in 1:length(points_f)){
  for(j in 1:length(points_f)){
    if(db8_pts_f$range[i] != db8_pts_f$range[j]){
      valley_matrix_f[i,j] <- (min(db8_pts_f$elev_ALOS[i], db8_pts_f$elev_ALOS[j]))
    }
  }
}
valley_barrier_f <- cbind(points_f, data.frame(valley_matrix_f))
names(valley_barrier_f)[1] <- "point"

# Pasture points
#   distance matrix for montane barriers
montane_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
for(i in 1:length(points_p)){
  for(j in 1:length(points_p)){
    if(db8_pts_p$east[i] != db8_pts_p$east[j]){
      montane_matrix_p[i,j] <- (4100 - max(db8_pts_p$elev_ALOS[i], db8_pts_p$elev_ALOS[j]))
    }
  }
}
montane_barrier_p <- cbind(points_p, data.frame(montane_matrix_p))
names(montane_barrier_p)[1] <- "point"

# distance matrix for valley barriers
valley_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
for(i in 1:length(points_p)){
  for(j in 1:length(points_p)){
    if(db8_pts_p$range[i] != db8_pts_p$range[j]){
      valley_matrix_p[i,j] <- (min(db8_pts_p$elev_ALOS[i], db8_pts_p$elev_ALOS[j]))
    }
  }
}
valley_barrier_p <- cbind(points_p, data.frame(valley_matrix_p))
names(valley_barrier_p)[1] <- "point"

gdmTab_f <- formatsitepair(db_format3_f, bioFormat=3, XColumn="lon", YColumn="lat",
                           siteColumn="point", 
                           predData = as.data.frame(db8_pts_f[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")]),
                           distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
gdmTab_p <- formatsitepair(db_format3_p, bioFormat=3, XColumn="lon", YColumn="lat",
                           sppColumn="scientificName", siteColumn="point", 
                           predData = as.data.frame(db8_pts_p[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")]),
                           distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))

# Fit GDM for forest points
#   Observed data
forest_gdm_obs <- gdm(gdmTab_f, geo = T)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_f <- gdmTab_f
forest_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
  pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
  withweights_f$weights <- pair_weights_f
  forest_gdm_bb[[i]] <- gdm(withweights_f, geo = T)
}

# Fit GDM for pasture points
#   Observed data
pasture_gdm_obs <- gdm(gdmTab_p, geo = T)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_p <- gdmTab_p
pasture_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
  pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
  withweights_p$weights <- pair_weights_p
  pasture_gdm_bb[[i]] <- gdm(withweights_p, geo = T)
}

gdms_raw_turnover <- list(forest_obs = forest_gdm_obs, pasture_obs = pasture_gdm_obs, forest_bb = forest_gdm_bb, pasture_bb = pasture_gdm_bb)
saveRDS(gdms_raw_turnover, file = "./gdm/gdms_raw_turnover.RDS")

######################### Modeled data: Turnover ###############################

species_predictions <- readRDS("./gdm/species_predictions_gdm.RDS")

a_rep_df <- bind_rows(species_predictions)

# Nombres de las columnas de abundancia
abun_cols <- grep("^abun__draw_", colnames(a_rep_df), value = TRUE)

# Crear lista con un dataframe por cada columna de abundancia
a_rep <- lapply(abun_cols, function(col_name) {
  df <- a_rep_df[, setdiff(colnames(a_rep_df), abun_cols)] # todas las columnas excepto abun__draw_XX
  df$a <- a_rep_df[[col_name]]  # nueva columna 'a' con la abundancia específica
  return(df)
})

# Asignar nombres a la lista
names(a_rep) <- abun_cols

forest_gdm_bb <- pasture_gdm_bb <- list()
for(k in 1:100){ # k indexes the posterior iteration
  print(k)
  db9 <- db7[a_rep[[k]]$a >= 1, ]
  db9_pts <- db9[!duplicated(db9$point), c("point", "lat", "lon", "elev_ALOS", "pasture", "precip","precip_ceccherini", "precip_IDEAM", "east", "range")]
  db9_f <- db9[db9$pasture == 0, c("scientificName", "point", "lat", "lon")]
  db9_pts_f <- db9_pts_f_stable <- db9_pts[db9_pts$pasture == 0, names(db9_pts) != "pasture"]
  db9_p <- db9[db9$pasture == 1, c("scientificName", "point", "lat", "lon")]
  db9_pts_p <- db9_pts_p_stable <- db9_pts[db9_pts$pasture == 1, names(db9_pts) != "pasture"]
  db9pts_f <- db9_pts_f_stable[db9_pts_f_stable$point %in% db9_f$point,] # exclude points with no detections
  db9pts_p <- db9_pts_p_stable[db9_pts_p_stable$point %in% db9_p$point,] # exclude points with no detections
  # forest
  points_f <- db9pts_f$point
  montane_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
  for(i in 1:length(points_f)){
    for(j in 1:length(points_f)){
      if(db9pts_f$east[i] != db9pts_f$east[j]){
        montane_matrix_f[i,j] <- (4100 - max(db9pts_f$elev_ALOS[i], db9pts_f$elev_ALOS[j]))
      }
    }
  }
  montane_barrier_f <- cbind(points_f, data.frame(montane_matrix_f))
  names(montane_barrier_f)[1] <- "point"
  
  valley_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
  for(i in 1:length(points_f)){
    for(j in 1:length(points_f)){
      if(db9pts_f$range[i] != db9pts_f$range[j]){
        valley_matrix_f[i,j] <- (min(db9pts_f$elev_ALOS[i], db9pts_f$elev_ALOS[j]))
      }
    }
  }
  valley_barrier_f <- cbind(points_f, data.frame(valley_matrix_f))
  names(valley_barrier_f)[1] <- "point"
  
  # pasture
  points_p <- db9pts_p$point
  montane_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
  for(i in 1:length(points_p)){
    for(j in 1:length(points_p)){
      if(db9pts_p$east[i] != db9pts_p$east[j]){
        montane_matrix_p[i,j] <- (4100 - max(db9pts_p$elev_ALOS[i], db9pts_p$elev_ALOS[j]))
      }
    }
  }
  montane_barrier_p <- cbind(points_p, data.frame(montane_matrix_p))
  names(montane_barrier_p)[1] <- "point"
  
  
  valley_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
  for(i in 1:length(points_p)){
    for(j in 1:length(points_p)){
      if(db9pts_p$range[i] != db9pts_p$range[j]){
        valley_matrix_p[i,j] <- (min(db9pts_p$elev_ALOS[i], db9pts_p$elev_ALOS[j]))
      }
    }
  }
  valley_barrier_p <- cbind(points_p, data.frame(valley_matrix_p))
  names(valley_barrier_p)[1] <- "point"
  
  # Format site-pair tables for forest and pasture
  db_spp_f <- unique(db9_f$scientificName)
  db_matrix_f <- matrix(data=0, nrow=nrow(db9_pts_f), ncol = length(db_spp_f))
  db_matrix_f[cbind(match(db9_f$point, db9_pts_f$point), match(db9_f$scientificName, db_spp_f))] <- 1
  forest_betapart <- beta.pair(db_matrix_f)
  db_format3_f <- cbind(db9_pts_f$point, as.matrix(forest_betapart$beta.sim))
  colnames(db_format3_f)[1] <- "point"
  db_format3_f[, 1] <- 1:nrow(db_format3_f)
  db_format3_f <- as.data.frame(db_format3_f)
  db_format3_f$point <- as.numeric(db_format3_f$point)
  
  db_spp_p <- unique(db9_p$scientificName)
  db_matrix_p <- matrix(data=0, nrow=nrow(db9_pts_p), ncol = length(db_spp_p))
  db_matrix_p[cbind(match(db9_p$point, db9_pts_p$point), match(db9_p$scientificName, db_spp_p))] <- 1
  pasture_betapart <- beta.pair(db_matrix_p)
  db_format3_p <- cbind(db9_pts_p$point, as.matrix(pasture_betapart$beta.sim))
  colnames(db_format3_p)[1] <- "point"
  db_format3_p[, "point"] <- 1:nrow(db_format3_p)
  db_format3_p <- as.data.frame(db_format3_p)
  db_format3_p$point <- as.numeric(db_format3_p$point)
  

  gdmTab_f <- formatsitepair(db_format3_f, bioFormat=3, XColumn="lon", YColumn="lat",
                             sppColumn="scientificName", siteColumn="point", 
                             predData = db9_pts_f[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                             distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
  gdmTab_p <- formatsitepair(db_format3_p, bioFormat=3, XColumn="lon", YColumn="lat",
                             sppColumn="scientificName", siteColumn="point", 
                             predData = db9_pts_p[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                             distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))
  
  # Fit GDM for forest points
  #   Observed data
  forest_gdm_obs <- gdm(gdmTab_f, geo = T)
  #   bayesian bootstrap resample
  n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
  n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
  s_all <- unique(c(n1,n2))
  s1 <- match(n1, s_all)
  s2 <- match(n2, s_all)
  withweights_f <- gdmTab_f
  
  
  site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
  pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
  withweights_f$weights <- pair_weights_f
  forest_gdm_bb[[k]] <- gdm(withweights_f, geo = T)
  
  # Fit GDM for pasture points
  #   Observed data
  pasture_gdm_obs <- gdm(gdmTab_p, geo = T)
  #   bayesian bootstrap resample
  n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
  n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
  s_all <- unique(c(n1,n2))
  s1 <- match(n1, s_all)
  s2 <- match(n2, s_all)
  withweights_p <- gdmTab_p
  
  
  site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
  pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
  withweights_p$weights <- pair_weights_p
  pasture_gdm_bb[[k]] <- gdm(withweights_p, geo = T)
  
}

gdms_modeled_simpsons <- list(forest_gdm_rep_bb = forest_gdm_bb, pasture_gdm_rep_bb = pasture_gdm_bb)
saveRDS(gdms_modeled_simpsons, file = "./gdm/gdms_modeled_simpsons_v6.RDS")

################################################################################
####################### Raw data gdm: Raup-crick ###############################
db8_pts <- readRDS("./gdm/db8_pts.RDS")
db8_pts_f <- db8_pts_f_stable <- db8_pts[db8_pts$pasture == 0, names(db8_pts) != "pasture"] # the "stable" version will be a version to come back to when we start modifying birdpts_f later on
db8_pts_p <- db8_pts_p_stable <- db8_pts[db8_pts$pasture == 1, names(db8_pts) != "pasture"]
db8_f <- db8[db8$pasture == 0, c("scientificName", "point", "lat", "lon")]
db8_p <- db8[db8$pasture == 1, c("scientificName", "point", "lat", "lon")]
db8_pts_p <- db8_pts_p[db8_pts_p$point %in% db8_p$point,] # exclude points with no detections in the raw pasture data (there are no such points in the forest data)
points_f <- db8_pts_f$point
points_p <- db8_pts_p$point

# Format site-pair tables for forest and pasture
db_spp_f <- unique(db8_f$scientificName)
db_matrix_f <- matrix(data=0, nrow=nrow(db8_pts_f), ncol = length(db_spp_f))
db_matrix_f[cbind(match(db8_f$point, db8_pts_f$point), match(db8_f$scientificName, db_spp_f))] <- 1
forest_raup <- vegan::raupcrick(db_matrix_f)
db_format3_f <- cbind(db8_pts_f$point, as.matrix(forest_raup))
colnames(db_format3_f)[1] <- "point"
db_format3_f[, 1] <- 1:nrow(db_format3_f)
db_format3_f <- as.data.frame(db_format3_f)
db_format3_f$point <- as.numeric(db_format3_f$point)

db_spp_p <- unique(db8_p$scientificName)
db_matrix_p <- matrix(data=0, nrow=nrow(db8_pts_p), ncol = length(db_spp_p))
db_matrix_p[cbind(match(db8_p$point, db8_pts_p$point), match(db8_p$scientificName, db_spp_p))] <- 1
pasture_raup <- vegan::raupcrick(db_matrix_p)
db_format3_p <- cbind(db8_pts_p$point, as.matrix(pasture_raup))
colnames(db_format3_p)[1] <- "point"
db_format3_p[, "point"] <- 1:nrow(db_format3_p)
db_format3_p <- as.data.frame(db_format3_p)
db_format3_p$point <- as.numeric(db_format3_p$point)

# For forest points:
#   distance matrix for montane barriers
montane_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
for(i in 1:length(points_f)){
  for(j in 1:length(points_f)){
    if(db8_pts_f$east[i] != db8_pts_f$east[j]){
      montane_matrix_f[i,j] <- (4100 - max(db8_pts_f$elev_ALOS[i], db8_pts_f$elev_ALOS[j]))
    }
  }
}
montane_barrier_f <- cbind(points_f, data.frame(montane_matrix_f))
names(montane_barrier_f)[1] <- "point"

#   distance matrix for valley barriers
valley_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
for(i in 1:length(points_f)){
  for(j in 1:length(points_f)){
    if(db8_pts_f$range[i] != db8_pts_f$range[j]){
      valley_matrix_f[i,j] <- (min(db8_pts_f$elev_ALOS[i], db8_pts_f$elev_ALOS[j]))
    }
  }
}
valley_barrier_f <- cbind(points_f, data.frame(valley_matrix_f))
names(valley_barrier_f)[1] <- "point"

# Pasture points
#   distance matrix for montane barriers
montane_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
for(i in 1:length(points_p)){
  for(j in 1:length(points_p)){
    if(db8_pts_p$east[i] != db8_pts_p$east[j]){
      montane_matrix_p[i,j] <- (4100 - max(db8_pts_p$elev_ALOS[i], db8_pts_p$elev_ALOS[j]))
    }
  }
}
montane_barrier_p <- cbind(points_p, data.frame(montane_matrix_p))
names(montane_barrier_p)[1] <- "point"

# distance matrix for valley barriers
valley_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
for(i in 1:length(points_p)){
  for(j in 1:length(points_p)){
    if(db8_pts_p$range[i] != db8_pts_p$range[j]){
      valley_matrix_p[i,j] <- (min(db8_pts_p$elev_ALOS[i], db8_pts_p$elev_ALOS[j]))
    }
  }
}
valley_barrier_p <- cbind(points_p, data.frame(valley_matrix_p))
names(valley_barrier_p)[1] <- "point"

gdmTab_f <- formatsitepair(db_format3_f, bioFormat=3, XColumn="lon", YColumn="lat",
                           siteColumn="point", 
                           predData = as.data.frame(db8_pts_f[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")]),
                           distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
gdmTab_p <- formatsitepair(db_format3_p, bioFormat=3, XColumn="lon", YColumn="lat",
                           sppColumn="scientificName", siteColumn="point", 
                           predData = as.data.frame(db8_pts_p[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")]),
                           distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))

# Fit GDM for forest points
#   Observed data
forest_gdm_obs <- gdm(gdmTab_f, geo = T)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_f <- gdmTab_f
forest_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
  pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
  withweights_f$weights <- pair_weights_f
  forest_gdm_bb[[i]] <- gdm(withweights_f, geo = T)
}

# Fit GDM for pasture points
#   Observed data
pasture_gdm_obs <- gdm(gdmTab_p, geo = T)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_p <- gdmTab_p
pasture_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
  pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
  withweights_p$weights <- pair_weights_p
  pasture_gdm_bb[[i]] <- gdm(withweights_p, geo = T)
}

gdms_raw_raup <- list(forest_obs = forest_gdm_obs, pasture_obs = pasture_gdm_obs, forest_bb = forest_gdm_bb, pasture_bb = pasture_gdm_bb)
saveRDS(gdms_raw_raup, file = "./gdm/gdms_raw_raup.RDS")

##### Raw data gdm: Raup-crick; no dist #####
# Fit GDM for forest points
#   Observed data
forest_gdm_obs <- gdm(gdmTab_f, geo = F)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_f <- gdmTab_f
forest_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
  pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
  withweights_f$weights <- pair_weights_f
  forest_gdm_bb[[i]] <- gdm(withweights_f, geo = F)
}

# Fit GDM for pasture points
#   Observed data
pasture_gdm_obs <- gdm(gdmTab_p, geo = F)
#   bayesian bootstrap resample
n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
s_all <- unique(c(n1,n2))
s1 <- match(n1, s_all)
s2 <- match(n2, s_all)
withweights_p <- gdmTab_p
pasture_gdm_bb <- list()
for(i in 1:100){
  print(i)
  site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
  pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
  withweights_p$weights <- pair_weights_p
  pasture_gdm_bb[[i]] <- gdm(withweights_p, geo = F)
}

gdms_raw_raup <- list(forest_obs = forest_gdm_obs, pasture_obs = pasture_gdm_obs, forest_bb = forest_gdm_bb, pasture_bb = pasture_gdm_bb)
saveRDS(gdms_raw_raup, file = "./gdm/gdms_raw_raup_nd.RDS")

######################### Modeled data: Raup-crick ###############################

species_predictions <- readRDS("./gdm/species_predictions_gdm.RDS")

a_rep_df <- bind_rows(species_predictions)

# Nombres de las columnas de abundancia
abun_cols <- grep("^abun__draw_", colnames(a_rep_df), value = TRUE)

# Crear lista con un dataframe por cada columna de abundancia
a_rep <- lapply(abun_cols, function(col_name) {
  df <- a_rep_df[, setdiff(colnames(a_rep_df), abun_cols)] # todas las columnas excepto abun__draw_XX
  df$a <- a_rep_df[[col_name]]  # nueva columna 'a' con la abundancia específica
  return(df)
})

# Asignar nombres a la lista
names(a_rep) <- abun_cols

forest_gdm_bb <- pasture_gdm_bb <- list()
for(k in 1:100){ # k indexes the posterior iteration
  print(k)
  db9 <- db7[a_rep[[k]]$a >= 1, ]
  db9_pts <- db9[!duplicated(db9$point), c("point", "lat", "lon", "elev_ALOS", "pasture", "precip","precip_ceccherini", "precip_IDEAM", "east", "range")]
  db9_f <- db9[db9$pasture == 0, c("scientificName", "point", "lat", "lon")]
  db9_pts_f <- db9_pts_f_stable <- db9_pts[db9_pts$pasture == 0, names(db9_pts) != "pasture"]
  db9_p <- db9[db9$pasture == 1, c("scientificName", "point", "lat", "lon")]
  db9_pts_p <- db9_pts_p_stable <- db9_pts[db9_pts$pasture == 1, names(db9_pts) != "pasture"]
  db9pts_f <- db9_pts_f_stable[db9_pts_f_stable$point %in% db9_f$point,] # exclude points with no detections
  db9pts_p <- db9_pts_p_stable[db9_pts_p_stable$point %in% db9_p$point,] # exclude points with no detections
  # forest
  points_f <- db9pts_f$point
  montane_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
  for(i in 1:length(points_f)){
    for(j in 1:length(points_f)){
      if(db9pts_f$east[i] != db9pts_f$east[j]){
        montane_matrix_f[i,j] <- (4100 - max(db9pts_f$elev_ALOS[i], db9pts_f$elev_ALOS[j]))
      }
    }
  }
  montane_barrier_f <- cbind(points_f, data.frame(montane_matrix_f))
  names(montane_barrier_f)[1] <- "point"
  
  valley_matrix_f <- matrix(data = 0, nrow = length(points_f), ncol = length(points_f))
  for(i in 1:length(points_f)){
    for(j in 1:length(points_f)){
      if(db9pts_f$range[i] != db9pts_f$range[j]){
        valley_matrix_f[i,j] <- (min(db9pts_f$elev_ALOS[i], db9pts_f$elev_ALOS[j]))
      }
    }
  }
  valley_barrier_f <- cbind(points_f, data.frame(valley_matrix_f))
  names(valley_barrier_f)[1] <- "point"
  
  # pasture
  points_p <- db9pts_p$point
  montane_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
  for(i in 1:length(points_p)){
    for(j in 1:length(points_p)){
      if(db9pts_p$east[i] != db9pts_p$east[j]){
        montane_matrix_p[i,j] <- (4100 - max(db9pts_p$elev_ALOS[i], db9pts_p$elev_ALOS[j]))
      }
    }
  }
  montane_barrier_p <- cbind(points_p, data.frame(montane_matrix_p))
  names(montane_barrier_p)[1] <- "point"
  
  
  valley_matrix_p <- matrix(data = 0, nrow = length(points_p), ncol = length(points_p))
  for(i in 1:length(points_p)){
    for(j in 1:length(points_p)){
      if(db9pts_p$range[i] != db9pts_p$range[j]){
        valley_matrix_p[i,j] <- (min(db9pts_p$elev_ALOS[i], db9pts_p$elev_ALOS[j]))
      }
    }
  }
  valley_barrier_p <- cbind(points_p, data.frame(valley_matrix_p))
  names(valley_barrier_p)[1] <- "point"
  
  # Format site-pair tables for forest and pasture
  db_spp_f <- unique(db9_f$scientificName)
  db_matrix_f <- matrix(data=0, nrow=nrow(db9_pts_f), ncol = length(db_spp_f))
  db_matrix_f[cbind(match(db9_f$point, db9_pts_f$point), match(db9_f$scientificName, db_spp_f))] <- 1
  forest_raup <- vegan::raupcrick(db_matrix_f)
  db_format3_f <- cbind(db9_pts_f$point, as.matrix(forest_raup))
  colnames(db_format3_f)[1] <- "point"
  db_format3_f[, 1] <- 1:nrow(db_format3_f)
  db_format3_f <- as.data.frame(db_format3_f)
  db_format3_f$point <- as.numeric(db_format3_f$point)
  
  db_spp_p <- unique(db9_p$scientificName)
  db_matrix_p <- matrix(data=0, nrow=nrow(db9_pts_p), ncol = length(db_spp_p))
  db_matrix_p[cbind(match(db9_p$point, db9_pts_p$point), match(db9_p$scientificName, db_spp_p))] <- 1
  pasture_raup <- vegan::raupcrick(db_matrix_p)
  db_format3_p <- cbind(db9_pts_p$point, as.matrix(pasture_raup))
  colnames(db_format3_p)[1] <- "point"
  db_format3_p[, "point"] <- 1:nrow(db_format3_p)
  db_format3_p <- as.data.frame(db_format3_p)
  db_format3_p$point <- as.numeric(db_format3_p$point)
  
  
  gdmTab_f <- formatsitepair(db_format3_f, bioFormat=3, XColumn="lon", YColumn="lat",
                             sppColumn="scientificName", siteColumn="point", 
                             predData = db9_pts_f[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                             distPreds = list(montane_barrier_f=montane_barrier_f, valley_barrier_f=valley_barrier_f))
  gdmTab_p <- formatsitepair(db_format3_p, bioFormat=3, XColumn="lon", YColumn="lat",
                             sppColumn="scientificName", siteColumn="point", 
                             predData = db9_pts_p[,c("point", "lat", "lon", "elev_ALOS", "precip_ceccherini")],
                             distPreds = list(montane_barrier_p=montane_barrier_p, valley_barrier_p=valley_barrier_p))
  
  # Fit GDM for forest points
  #   Observed data
  forest_gdm_obs <- gdm(gdmTab_f, geo = T)
  #   bayesian bootstrap resample
  n1 <- paste0(gdmTab_f$s1.xCoord, "_", gdmTab_f$s1.yCoord)
  n2 <- paste0(gdmTab_f$s2.xCoord, "_", gdmTab_f$s2.yCoord)
  s_all <- unique(c(n1,n2))
  s1 <- match(n1, s_all)
  s2 <- match(n2, s_all)
  withweights_f <- gdmTab_f
  
  
  site_weights_f <- rdirichlet(1, rep(1, length(points_f)))
  pair_weights_f <- site_weights_f[s1] * site_weights_f[s2]
  withweights_f$weights <- pair_weights_f
  forest_gdm_bb[[k]] <- gdm(withweights_f, geo = T)
  
  # Fit GDM for pasture points
  #   Observed data
  pasture_gdm_obs <- gdm(gdmTab_p, geo = T)
  #   bayesian bootstrap resample
  n1 <- paste0(gdmTab_p$s1.xCoord, "_", gdmTab_p$s1.yCoord)
  n2 <- paste0(gdmTab_p$s2.xCoord, "_", gdmTab_p$s2.yCoord)
  s_all <- unique(c(n1,n2))
  s1 <- match(n1, s_all)
  s2 <- match(n2, s_all)
  withweights_p <- gdmTab_p
  
  
  site_weights_p <- rdirichlet(1, rep(1, length(points_p)))
  pair_weights_p <- site_weights_p[s1] * site_weights_p[s2]
  withweights_p$weights <- pair_weights_p
  pasture_gdm_bb[[k]] <- gdm(withweights_p, geo = T)
  
}

gdms_modeled_raup <- list(forest_gdm_rep_bb = forest_gdm_bb, pasture_gdm_rep_bb = pasture_gdm_bb)
saveRDS(gdms_modeled_raup, file = "./gdm/gdms_modeled_raup_v6.RDS")

################################################################################
library(readxl)
library(reticulate)

db1 <- read_excel("abundance/Scarabaeinae_database_2024.xlsx", sheet = "Scarabaeinae_database_2024")

unique(db1$genus)
length(unique(db1$genus))
length(unique(db7$point))
length(unique(db7$cluster))

##### Set up GEE session #####
reticulate::use_condaenv("r-reticulate", required = TRUE)
library(rgee)
ee$Initialize()
np <- import("numpy")
pd <- import("pandas")


# Load raster files
ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
SRTM30 <- ee$Image("USGS/SRTMGL1_003")
RESOLVE <- ee$FeatureCollection("RESOLVE/ECOREGIONS/2017")
reg_ras <- RESOLVE$reduceToImage(reducer = ee$Reducer$first(),properties = list("ECO_ID"))$unmask(0)$reproject('epsg:4326',NULL,ee$Number(30.922080775909325))
biom_ras <- RESOLVE$reduceToImage(reducer = ee$Reducer$first(),properties = list("BIOME_NUM"))$unmask(0)$reproject('epsg:4326',NULL,ee$Number(30.922080775909325))

##### Extract raster values from all points #####
# Featurecollection of point geometries
pts_ALOS <- ALOS$select('AVE_DSM')$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
ALOSelev <- sapply(c(1:length(pts_ALOS$features)),function(x)pts_ALOS$features[[x]]$properties$mean)

pts_SRTM <- SRTM30$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
SRTMelev <- sapply(c(1:length(pts_SRTM$features)),function(x)pts_SRTM$features[[x]]$properties$mean)

geompts <- sapply(1:nrow(pts),function(x)ee$Geometry$Point(c(pts$lon[x],pts$lat[x])))
geompts <- ee$FeatureCollection(c(unlist(geompts)))

pts_REG <- reg_ras$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
ecoreg <- sapply(c(1:length(pts_REG$features)),function(x)pts_REG$features[[x]]$properties$mean)

pts_BIOM <- biom_ras$reduceRegions(geompts, ee$Reducer$mean())$getInfo()
biome <- sapply(c(1:length(pts_BIOM$features)),function(x)pts_BIOM$features[[x]]$properties$mean)

# Combined dataframe
spatialdata <- cbind.data.frame(point_id=pts$point_id,ALOSelev,SRTMelev,ecoreg,biome)

################################################################################
library(rnaturalearth)
library(ggplot2)
study_area <- st_read("F:/Capas/America/ecoregions/ecoreg.shp")
study_area <- st_make_valid(study_area)

db7_pts <- db7[!duplicated(db7$point), c("subregion", "cluster","point", "lat", "lon", "elev_ALOS", "pasture")]

db7_pts_sf <- st_as_sf(db7_pts, coords = c("lon", "lat"), crs = 4326)

# Realizar el join espacial
db7_joined <- st_join(db7_pts_sf, study_area[, c("ECO_NAME", "BIOME_NAME")])

# Convertir de vuelta a data.frame si quieres
db7_final <- as.data.frame(db7_joined)

colombia <- ne_countries(scale = "medium", country = "Colombia", returnclass = "sf")

ggplot() +
  geom_sf(data = colombia, fill = "gray95", color = "black") +  # fondo de Colombia
  geom_sf(data = study_area, aes(fill = ECO_NAME), alpha = 0.4, color = NA) +  # ecorregiones
  geom_sf(data = db7_pts_sf, color = "red", size = 2) +  # puntos
  theme_bw() +
  labs(title = "Puntos de muestreo y ecorregiones en Colombia",
       fill = "Ecorregión")


db7_final <- db7_final %>%
  mutate(
    ECO_NAME = ifelse(point %in% c(paste0("SGP", 4:9), paste0("PSP", 1:7)), 
                      "Caqueta moist forests", 
                      ifelse(point == "ORP2", 
                             "Cordillera Oriental montane forests", 
                             ifelse(subregion == "subregion_llanosWest",
                                    "Apure-Villavicencio Dry Forests",
                                    ifelse(subregion == "subregion_paujil",
                                           "Magdalena-Urabá Moist Forests", ECO_NAME)))),
    
    BIOME_NAME = ifelse(point %in% c(paste0("SGP", 4:9), paste0("PSP", 1:7)), 
                        "Tropical & Subtropical Moist Broadleaf Forests", 
                        ifelse(point == "ORP2",
                               "Tropical & Subtropical Moist Broadleaf Forests", 
                               ifelse(subregion == "subregion_llanosWest",
                                      "Tropical and Subtropical Dry Broadleaf Forests",
                                      ifelse(subregion == "subregion_paujil",
                                             "Tropical and Subtropical Moist Broadleaf Forests", BIOME_NAME))))
  )


eco_point <- db7_final %>%
  st_drop_geometry() %>%
  group_by(ECO_NAME) %>%
  summarise(
    n_clusters = n_distinct(cluster),
    n_puntos = n(),
    .groups = "drop")
