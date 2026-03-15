# ------------------------------------------------------------------------------
# Spatial prediction of dung beetle abundance across Colombia
#
# This script generates spatial predictions of species-specific dung beetle
# abundance across mainland Colombia using the Bayesian hierarchical model
# fitted in module `7_models`.
#
# Predictions are produced on a 2 × 2 km grid covering mainland Colombia.
# Elevation data are retrieved from Google Earth Engine (ALOS World 3D DEM),
# and grid cells are grouped into spatial subregions matching the spatial
# scale used in field sampling.
#
# Species-level predictors, including functional traits and elevation
# preferences, are integrated to construct the prediction dataset. Two
# alternative land-use scenarios are simulated for each grid cell:
# forest (pasture = 0) and pasture (pasture = 1).
#
# Abundance predictions are obtained from posterior draws of the Bayesian
# hierarchical negative binomial model implemented in `brms`. Posterior
# linear predictors are calculated for each species and grid cell and used
# to simulate abundance values while propagating model uncertainty.
# ------------------------------------------------------------------------------

library(sf)
library(rgee)
library(dplyr)
library(brms)
library(posterior)
setwd("C:/Users/PC/Dropbox/CO_DBdata")

######################## Spatial prediction grid ###############################
# Define projection and extract elevation values across mainland Colombia
# using Google Earth Engine (ALOS World 3D Digital Surface Model).
#
# Elevation is summarized for a 2 km grid covering the country.

`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

###################### Extract elevation from Google Earth Engine ##############
# Retrieve elevation values and geographic coordinates for grid cells
# across Colombia.
reticulate::use_condaenv("r-reticulate", required = TRUE)
# ee_Authenticate()
ee_Initialize()

roi <- ee$FeatureCollection('USDOS/LSIB_SIMPLE/2017')$filterMetadata('country_na', 'equals', 'Colombia')

DEM <- ee$Image("JAXA/ALOS/AW3D30/V2_2")$select(list("AVE_DSM"),list("elevation"))

latlng <- ee$Image$pixelLonLat()$addBands(DEM)$reduceRegion(reducer = ee$Reducer$toList(),
                                                            geometry = roi,
                                                            maxPixels = 10^9,
                                                            scale=2000)

lats <- ee$Array(latlng$get("latitude"))$getInfo()
lngs <- ee$Array(latlng$get("longitude"))$getInfo()
elevs <- ee$Array(latlng$get("elevation"))$getInfo()

######################## Convert to spatial dataset ############################
# Convert extracted coordinates and elevation values into a spatial dataset
# used as the basis for prediction points.
elevation1 <- data.frame(lon = lngs, lat = lats, elev = elevs) |>
  st_as_sf(coords = c("lon","lat"), remove = FALSE) |>
  st_set_crs(4326)

######################## Define spatial subregions #############################
# Group grid cells into spatial subregions using latitude–longitude blocks.
# These subregions match the spatial scale used in field sampling and
# correspond to the hierarchical structure included in the Bayesian model.
# add subregion assignments

u_lats <- unique(lats)[order(unique(lats))]
u_lngs <- unique(lngs)[order(unique(lngs))]

get_lat_block <- Vectorize(function(lat) {
  floor(max(which(u_lats <= lat)) / 10)
})
get_lon_block <- Vectorize(function(lon) {
  floor(max(which(u_lngs <= lon)) / 10)
})

elevation <- elevation1 |>
  mutate(
    lat_block = get_lat_block(lats),
    lon_block = get_lon_block(lngs),
  ) |>
  mutate(
    subregion = paste0("sr_", lat_block, "__", lon_block)
    )

######################## Species geographic filtering ##########################
# Restrict prediction grid cells to locations within each species'
# geographic range estimated in module `3_Geographic_range`.
db_range <- readRDS("./geographic_range/geographic_range.rds")

species_frames <- list()

for(i in 1:nrow(db_range)) {
  print(i)
  species_frames[[i]] <-
    elevation[db_range[i,],] |>
    mutate(scientificName = db_range$scientific[i])
}

scale_elev <- function(e, l, u) {
  2 * (e - l) / (u - l) - 1
}

######################## Integrate species traits ##############################
# Add species-level ecological traits and elevation limits used to calculate
# standardized elevation predictors for the model.

db_traits <- readRDS("./abundance/db5_traits.RDS") |>
  filter(!duplicated(scientificName)) |>
  select(scientificName, sp_elev_lower2, sp_elev_upper2, sp_region_amazon,
         sp_region_llanos, sp_region_caribbean, sp_region_snsm, sp_region_andes,
         sp_region_eastern, sp_slope_ECe, sp_slope_ECw, sp_slope_CCe, sp_slope_CCw,
         sp_slope_WCe, sp_slope_WCw, sp_slope_SNSM, nest_guild, diet_range, activity,
         bodysize, legratio)

pd <- do.call(rbind, species_frames) |>
  full_join(db_traits, by = "scientificName") |>
  mutate(
    elev_standard = scale_elev(elev, sp_elev_lower2, sp_elev_upper2),
    elev_standard_squared = elev_standard^2
    )

pd_pasture <- pd |>
  mutate(pasture = 1)
pd_forest <- pd |>
  mutate(pasture = 0)

######################## Build prediction dataset ##############################
# Combine spatial grid cells with species traits and environmental predictors.
# Two prediction scenarios are generated:
#   • forest (pasture = 0)
#   • pasture (pasture = 1)

prediction_data <- rbind(pd_pasture, pd_forest) |>
  mutate(subregion_species = paste0(subregion, "__", scientificName))

saveRDS(prediction_data, "./prediction_data.RDS")
prediction_data <- readRDS("prediction_data.rds")

######################## Load Bayesian hierarchical model ######################
# Import the fitted Bayesian abundance model generated in module `7_models`.
mod <- readRDS("./db_mod_abundance.rds")
mod_draws <- as_draws_df(mod)

######################## Posterior linear predictions ##########################
# Compute posterior linear predictors for each species and grid cell
# using sampled posterior draws from the Bayesian model.
species_preds <- cluster_sds <- phis <- list()
draw_ids <- 30 * c(1:100)
for(i in 1:nrow(db_traits)){
  print(i)
  species_preds[[i]] <- posterior_linpred(
    mod, 
    draw_ids = draw_ids,
    re_formula = ~ (1 + pasture + elev_standard + elev_standard_squared | scientificName) +
      (1 | subregion_species),
    newdata = prediction_data |>
      filter(scientificName == db_traits$scientificName[i]),
    allow_new_levels = TRUE,
    sample_new_levels = "gaussian"
  )
  cluster_sds[[i]] <- mod_draws$sd_cluster_species__Intercept[draw_ids]
  phis[[i]] <- mod_draws$shape[draw_ids]
}

######################## Abundance simulation ##################################
# Simulate abundance values from the negative binomial distribution while
# incorporating posterior uncertainty and cluster-level variability.
simulate_abun <- Vectorize(function(linpred, cluster_sd, phi) {
  out <- 0
  for(i in 1:16){
    lp <- linpred + rnorm(1, 0, cluster_sd)
    abun <- rnbinom(5, size = phi, mu = exp(lp))
    out <- out + mean(abun)/16
  }
  out
})

######################## Species-level abundance predictions ###################
# Generate predicted abundances for each species across the spatial grid
# under both land-use scenarios.

species_predictions <- list()
for(i in 1:nrow(db_traits)){
  print(i)
  pd_species <- prediction_data |>
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

######################## Save prediction outputs ###############################
# Save simulated abundance predictions for all species.

saveRDS(species_predictions, "species_predictions_100draws.RDS")
