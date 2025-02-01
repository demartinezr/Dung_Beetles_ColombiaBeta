library(sf)
library(dplyr)


setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata")

db_predictions <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/species_predictions.rds")
df_predictions <- bind_rows(lapply(db_predictions, function(df) {
  df %>% select(scientificName, pasture, starts_with("abun__draw_"), geometry)
}))

df_predictions <- df_predictions %>%
      rowwise() %>%
      mutate(median_abundance = median(c_across(starts_with("abun__draw_")), na.rm = TRUE),
      mean_abundance = mean(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
      ungroup()



saveRDS(df_predictions, "df_predictions.rds")
gc()
df_predictions <- sf_combined %>%
  select(scientificName, pasture, starts_with("abun__draw_"), geometry)
gc()

df_predictions <- df_predictions %>%
  rowwise() %>%
  mutate(median_abundance = median(c_across(starts_with("abun__draw_")), na.rm = TRUE),
         mean_abundance = mean(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
  ungroup()

