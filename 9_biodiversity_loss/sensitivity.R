library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata")

db_predictions <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/species_predictions.rds")
################################################################################
## get multiplicative change of abudance by draw for each species
calculate_multiplicative_change <- function(df) {
  forest <- df %>% filter(pasture == 1)
  pasture <- df %>% filter(pasture == 0)
  
  abundance_ratio <- forest %>%
    select(starts_with("abun__draw_")) / pasture %>%
    select(starts_with("abun__draw_"))
  
  colnames(abundance_ratio) <- gsub("abun__draw_", "ratio__draw_", colnames(abundance_ratio))
  
  result_df <- forest %>%
    select(scientificName) %>%
    bind_cols(abundance_ratio)
  
  return(result_df)
}
ratio_draw <- lapply(db_predictions, calculate_multiplicative_change)

saveRDS(ratio_draw, "./ratio_draw.rds")
#
###############################################################################
# Global sensitivities and mean ratio by draw for each species (columns)
ratio_draw <- readRDS("./ratio_draw.rds")

mean_ratio <- function(df) {
  df %>%
    summarise(across(starts_with("ratio__draw_"), 
                     ~mean(.x[is.finite(.x)], na.rm = TRUE))) %>% 
    mutate(scientificName = df$scientificName[1])
}

mean_ratio_draw <- bind_rows(lapply(ratio_draw, mean_ratio))


global_ratio <- mean_ratio_draw %>%
  mutate(mean_ratio_draw = rowMeans(select(., starts_with("ratio__draw_")), na.rm = TRUE))

mean_global_ratio <- global_ratio %>%
  select(scientificName, mean_ratio_draw)

# mean_global_ratio <- mean_global_ratio %>% # ratio > 1000
#  filter(!scientificName %in% c("Sulcophanaeus_leander", "Coprophanaeus_edmondsi", "Deltochilum_orbiculare"))

ggplot(mean_global_ratio, aes(x = log10(mean_ratio_draw))) +
  geom_histogram(binwidth = 0.1, fill = "grey90", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +  # Línea en x = 1
  labs(
    x = "Sensitivity",
    y = "Frecuency",
    title = "The distribution of species-specific sensitivities to forest-pasture conversion"
  ) +
  theme_classic()

write.csv(mean_global_ratio, "./mean_global_ratio.csv")

###############################################################################
# Global sensitivities and mean ratio by species (rows)
ratio_draw <- readRDS("./ratio_draw.rds")
plan(multicore, workers = 4)
sensitivity <- future_lapply(ratio_draw, function(df) {
  selected_cols <- df %>% select(starts_with("ratio__draw_"))
  
  df %>%
    mutate(
      median_ratio = apply(selected_cols, 1, function(x) {
        x <- x[is.finite(x)]
        if (length(x) > 0) median(x, na.rm = TRUE) else NA
      }),
      mean_ratio = apply(selected_cols, 1, function(x) {
        x <- x[is.finite(x)]
        if (length(x) > 0) mean(x, na.rm = TRUE) else NA
      })
    )
}, future.packages = c("dplyr"), future.seed = TRUE)

saveRDS(sensitivity, "./sensitivity.rds")

# median and mean ratio by species
sensitivity <- readRDS("./sensitivity.rds")

sensitivity_col <- future_lapply(sensitivity, function(df) {
  mean_ratio <- df$mean_ratio[is.finite(df$mean_ratio)]
  median_ratio <- df$median_ratio[is.finite(df$median_ratio)]
  
  ratio_col_mean <- mean(mean_ratio, na.rm = TRUE)
  ratio_col_median <- median(median_ratio, na.rm = TRUE)
  
  result_df <- data.frame(
    scientificName = df$scientificName[1],
    ratio_col_mean = ratio_col_mean,
    ratio_col_median = ratio_col_median
  )
  
  return(result_df)
}) %>%
  bind_rows()  
#
## The distribution of species-specific sensitivities to forest conversion
sensitivity_col <- sensitivity_col %>%
  filter(!scientificName %in% c("Sulcophanaeus_leander", "Coprophanaeus_edmondsi"))


ggplot(sensitivity_col, aes(x = log10(ratio_col_mean))) +
  geom_histogram(binwidth = 0.1, fill = "gray90", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +  # Línea en x = 1
  labs(
    x = "sensitivity",
    y = "Frecuency",
    title = "The distribution of species-specific sensitivities to forest-pasture conversion"
  ) +
  theme_classic()
write.csv(sensitivity_col, "./sensitivity_col_row.csv")

###############################################################################
#
# mean and median abundance by species (rows) and global sensitivities 
#
db_predictions <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/species_predictions.rds")

# median and mean ratio for each species based on 10 iterations
db_sensitivity <- lapply(db_predictions, function(df) {
  df %>%
    rowwise() %>%
    mutate(
      median_abundance = median(c_across(starts_with("abun__draw_")), na.rm = TRUE),
      mean_abundance = mean(c_across(starts_with("abun__draw_")), na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(scientificName, pasture, median_abundance, mean_abundance, geometry)  
})
saveRDS(db_sensitivity, "./sp_sensitivity.rds")

# get mean and median multiplicative change of abundance for each species 
# (sensitivity N forest / N pasture)
db_sensitivity <- readRDS("./sp_sensitivity.rds")

plan(multicore, workers = 4)
sensitivity_sp <- future_lapply(db_sensitivity, function(df) {
  # pasture (0) and forest (1) scenarios
  pasture <- df %>% filter(pasture == 0)
  forest <- df %>% filter(pasture == 1)
  # Perform the ratio calculation
  mean_ratio <- forest$mean_abundance / pasture$mean_abundance
  median_ratio <- forest$median_abundance / pasture$median_abundance
  result_df <- forest %>%
    select(scientificName, geometry) %>%
    mutate(ratio_mean = mean_ratio, ratio_median = median_ratio)
  return(result_df)
})
# statistics for mean and median ratios
sensitivity_col <- future_lapply(sensitivity_sp, function(df) {
  ratio_mean <- df$ratio_mean[is.finite(df$ratio_mean)]
  ratio_median <- df$ratio_median[is.finite(df$ratio_median)]
  
  ratio_col_mean <- mean(ratio_mean, na.rm = TRUE)
  ratio_col_median <- median(ratio_median, na.rm = TRUE)
  
  result_df <- data.frame(
    scientificName = df$scientificName[1],
    ratio_col_mean = ratio_col_mean,
    ratio_col_median = ratio_col_median
  )
  
  return(result_df)
}) %>%
  bind_rows()  

sensitivity_col <- sensitivity_col %>%
  filter(!scientificName %in% c("Sulcophanaeus_leander", "Coprophanaeus_edmondsi"))

ggplot(sensitivity_col, aes(x = log10(ratio_col_mean))) +
  geom_histogram(binwidth = 0.1, fill = "grey90", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +  # Línea en x = 1
  labs(
    x = "sensitivity",
    y = "Frecuency",
    title = "The distribution of species-specific sensitivities to forest conversion"
  ) +
  theme_classic()
write.csv(sensitivity_col, "./sensitivity_col.csv")

###########################################################################
db_predictions <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/species_predictions.rds")

# Get the region names from the list names
sp_names <- names(db_predictions)

# Function to calculate the proportional change in abundance by species and region
abundance_change <- function(sf_df, sp_name) {
  sf_df %>%
    st_drop_geometry() %>% 
    # Divide into forest (pasture = 1) and grassland (pasture = 0)
    group_by(pasture) %>%
    summarise(across(starts_with("abun__draw_"), \(x) sum(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_wider(names_from = pasture, values_from = starts_with("abun__draw_"), 
                names_glue = "{.value}_pasture{pasture}") %>%
    mutate(across(ends_with("_pasture1"), 
                  ~ .x / get(sub("_pasture1", "_pasture0", cur_column())), 
                  .names = "ratio_{.col}")) %>%
    select(starts_with("ratio_")) %>%
    mutate(scientificName = sp_name)
}

# Apply the function to each element of the list
ratio_draw <- map2(db_predictions, sp_names, abundance_change)

mean_ratio <- function(df, df_name) {
  df %>%
    rowwise() %>% 
    mutate(
      mean_ratio = mean(c_across(starts_with("ratio_abun__draw_"))[is.finite(c_across(starts_with("ratio_abun__draw_")))], na.rm = TRUE),
      median_ratio = median(c_across(starts_with("ratio_abun__draw_"))[is.finite(c_across(starts_with("ratio_abun__draw_")))], na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(scientificName, mean_ratio, median_ratio) %>%
    mutate(scientificName = df_name)
}
mean_ratio_draw <- bind_rows(mapply(mean_ratio, ratio_draw, names(ratio_draw), SIMPLIFY = FALSE))
mean_ratio_draw1 <- mean_ratio_draw %>%
  filter(!scientificName %in% c("Sulcophanaeus_leander", "Coprophanaeus_edmondsi", "Deltochilum_orbiculare"))

ggplot(mean_ratio_draw1, aes(x = log10(mean_ratio))) +
  geom_histogram(binwidth = 0.1, fill = "grey90", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 1) +  # Línea en x = 1
  labs(
    x = "sensitivity",
    y = "Frecuency",
    title = "The distribution of species-specific sensitivities to forest conversion"
  ) +
  theme_classic()
