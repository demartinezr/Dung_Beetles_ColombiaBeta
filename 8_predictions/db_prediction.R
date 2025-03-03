library(sf)
library(rgee)
library(dplyr)
library(brms)
library(posterior)
setwd("C:/Users/PC/Dropbox/CO_DBdata")

##### extract elevations on a 2 km grid over mainland colombia #####
`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

##### Get elevation raster across Colombia #####
ee_Authenticate()
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

# Convert to spatial points
elevation1 <- data.frame(lon = lngs, lat = lats, elev = elevs) |>
  st_as_sf(coords = c("lon","lat"), remove = FALSE) |>
  st_set_crs(4326)

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

##### for each species subset dataframe and join #####
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
  mutate(pasture = 0)
pd_forest <- pd |>
  mutate(pasture = 1)

prediction_data <- rbind(pd_pasture, pd_forest) |>
  mutate(subregion_species = paste0(subregion, "__", scientificName))

saveRDS(prediction_data, "./prediction_data.RDS")
prediction_data <- readRDS("prediction_data.rds")

##### do prediction #####
mod <- readRDS("./db_mod_abundance.rds")
mod_draws <- as_draws_df(mod)

species_preds <- cluster_sds <- phis <- list()
draw_ids <- 50 * c(1:60)
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

saveRDS(species_predictions, "species_predictions_60draws.RDS")
