library(sf)
library(dplyr)
library(brms)
library(posterior)

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

# new data set for predictions
db_traits <- readRDS("./abundance/db5_traits.rds") 
prediction_data <- readRDS("prediction_data.rds")

##### do prediction #####
mod <- readRDS("./db_mod_abundance.rds")
mod_draws <- as_draws_df(mod)

species_preds <- cluster_sds <- phis <- list()
 draw_ids <- 100 * c(1:30)
# draw_ids <- 30 * c(1:100)
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

saveRDS(species_predictions, "species_predictions_30draws.rds")
