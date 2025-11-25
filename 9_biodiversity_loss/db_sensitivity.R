setwd("C:/Users/PC/Dropbox/CO_DBdata")

library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

#db_predictions <- readRDS("./species_predictions.rds")
db_predictions <- readRDS("./species_predictions_100.rds")
#ecoregions_all <- as.data.frame(readRDS("./ecoregions_100.rds"))

############# Analysis of species specific sensitivities #####################

# Get the region names from the list names
sp_names <- names(db_predictions)
# calculate the proportional change in abundance for species and region by draw
abundance_change <- function(sf_df, sp_name) {
  sf_df %>%
    st_drop_geometry() %>% 
    # Divide into forest (pasture = 0) and pasture (pasture = 1)
    group_by(pasture) %>%
    summarise(across(starts_with("abun__draw_"), \(x) sum(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_wider(names_from = pasture, values_from = starts_with("abun__draw_"), 
                names_glue = "{.value}_pasture{pasture}") %>%
    mutate(across(ends_with("_pasture1"), 
                  ~ get(sub("_pasture1", "_pasture0", cur_column())) / .x, 
                  .names = "ratio_{.col}")) %>%
    dplyr::select(starts_with("ratio_")) %>%
    mutate(scientificName = sp_name)
}
# Apply the function for each species
ratio_draw <- map2(db_predictions, sp_names, abundance_change)

# get mean and median for species
mean_ratio <- function(df, df_name) {
  df %>%
    rowwise() %>% 
    mutate(
      mean_ratio = mean(c_across(starts_with("ratio_abun__draw_"))[is.finite(c_across(starts_with("ratio_abun__draw_")))], na.rm = TRUE),
      median_ratio = median(c_across(starts_with("ratio_abun__draw_"))[is.finite(c_across(starts_with("ratio_abun__draw_")))], na.rm = TRUE)
    ) %>%
    ungroup() %>%
    dplyr::select(scientificName, mean_ratio, median_ratio) %>%
    mutate(scientificName = df_name)
}
# Apply the function to each species
mean_ratio_draw <- bind_rows(mapply(mean_ratio, ratio_draw, names(ratio_draw), SIMPLIFY = FALSE))

# plot for the distribution of the mean abundance change across species 
fig_1c <- ggplot(mean_ratio_draw, aes(x = log10(median_ratio))) +
  geom_histogram(binwidth = 0.2, fill = "grey90", color = "black", alpha = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_vline(xintercept = c(-1, 1, 2, 3, 4), 
             linetype = "dotted", color = "black", linewidth = 0.8) +
  
  scale_x_continuous(name = "Sensitivity",
    breaks = c(-1, 0, 1, 2, 3, 4),
    labels = c("10x pasture", "0", "10x forest", "100x forest", "1,000x forest", "10,000x forest")) +
  labs( 
    y = "Frequency") +
#    title = "The distribution of species-specific sensitivities to forest-pasture conversion") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 12, face = "bold"))

ggsave("./fig_1c.jpeg", plot = fig_1c, width = 8.5, height = 3, units = "in",        
       dpi = 300, device = "jpeg")
################# Analysis of sensitivities for ecoregions scale ####################
setwd("C:/Users/PC/Dropbox/CO_DBdata")
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggridges)
library(gridExtra)

# Abundance predictions for 243 species of dung beetles with 10 iterations in 
# pasture and forest
db_predictions <- readRDS("./species_predictions.rds")

# Study area ecoregions
study_area <- st_read("F:/Capas/America/ecoregions/ecoreg.shp")
study_area <- st_make_valid(study_area)
st_is_valid(study_area)
# Create a list to store sf dataframes by ecoregion
ec_points <- vector("list", length = nrow(study_area))
names(ec_points) <- study_area$ECO_NAME  # Use ecoregion names

# Iterate over each ecoregion
for (i in seq_len(nrow(study_area))) {
  cat("ecoregion:", i, "of", nrow(study_area), "\n") 
  polygon <- study_area[i, ]
  # Filter points of all species that fall within the ecoregion
  points_in_ecoregion <- do.call(rbind, lapply(db_predictions, function(df) {
    st_filter(df, polygon)
  }))
  # Store the result in the list
  ec_points[[i]] <- points_in_ecoregion
}
saveRDS(ec_points, "./Analysis/mean_abundance/ecoregions_predictions_2.rds")

########################## ecoregions analysis ###############################
ecoregions_predictions <- readRDS("./Analysis/mean_abundance/ecoregions_predictions_2.rds")

regional_richness <- data.frame(
  sf_dataframe = names(ecoregions_predictions), 
  Total = sapply(ecoregions_predictions, function(x) n_distinct(x$scientificName)), row.names = NULL 
)
# Get the region names from the list names
region_names <- gsub(" forests", "", names(ecoregions_predictions))

# Function to calculate the proportional change in species abundance by region
abundance_change <- function(ecoregions_predictions, region_name) {
  ecoregions_predictions %>%
    st_drop_geometry() %>% 
    # Divide into forest (pasture = 1) and pasture (pasture = 0)
    group_by(scientificName, pasture) %>%
    summarise(across(starts_with("abun__draw_"), \(x) sum(x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_wider(names_from = pasture, values_from = starts_with("abun__draw_"), 
                names_glue = "{.value}_pasture{pasture}") %>%
    mutate(across(ends_with("_pasture1"), 
                  ~ .x / get(sub("_pasture1", "_pasture0", cur_column())), 
                  .names = "ratio_{.col}")) %>%
    select(scientificName, starts_with("ratio_")) %>%
    mutate(region = region_name)
}

# Apply the function to each element of the list
ratio_draw <- map2(ecoregions_predictions, region_names, abundance_change)

mean_ratio <- function(df, df_name) {
  df %>%
    summarise(across(starts_with("ratio_abun__draw_"), 
                     list(mean = ~mean(.x[is.finite(.x)], na.rm = TRUE),
                          p25  = ~quantile(.x[is.finite(.x)], probs = 0.25, na.rm = TRUE),
                          p50  = ~quantile(.x[is.finite(.x)], probs = 0.50, na.rm = TRUE),
                          p75  = ~quantile(.x[is.finite(.x)], probs = 0.75, na.rm = TRUE)))) %>% 
    pivot_longer(cols = starts_with("ratio_abun__draw_"),
                 names_to = c("draw", ".value"),
                 names_pattern = "ratio_abun__draw_(\\d+)_(.*)") %>%
    mutate(ecoregion = df_name)
}
mean_ratio_draw <- bind_rows(mapply(mean_ratio, ratio_draw, names(ratio_draw), SIMPLIFY = FALSE))
mean_ratio_col <- readRDS("./mean_ratio_col_2.rds")
mean_ratio_draw$ecoregion <- gsub(" forests", "", mean_ratio_draw$ecoregion)
mean_ratio_draw <- bind_rows(mean_ratio_draw, mean_ratio_col)

plot_p25 <- ggplot(mean_ratio_draw, aes(x = pasture1_p25, y = reorder(ecoregion, -pasture1_p25), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.6, scale = 1.5) +  
  labs(title = "25th percentile",
       x = "sensitivity\n(N forest/N pasture)",
       y = "Ecoregion") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=6),
        axis.title.x = element_blank(),  
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  scale_x_log10(breaks = c(0.1, 0.5, 1, 3, 6), 
                labels = c("0.1","0.5", "1", "3", "6"),
                limits = c(0.1, 6)) +
  stat_density(geom = "line", position = "identity", aes(x = pasture1_p25), adjust = 1, kernel = "gaussian")

plot_mean <- ggplot(mean_ratio_draw, aes(x = pasture1_mean, y = reorder(ecoregion, -pasture1_mean), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.6, scale = 1.5) +  
  labs(title = "Mean of the multiplicative abundance change across species",
       x = "sensitivity\n(N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0, 50)) +
  stat_density(geom = "line", position = "identity", aes(x = pasture1_mean), adjust = 1, kernel = "gaussian")

plot_p50 <- ggplot(mean_ratio_draw, aes(x = pasture1_p50, y = reorder(ecoregion, -pasture1_p50), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.6, scale = 1.5) +  
  labs(title = "Median",
       x = "sensitivity\n(N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=6),
        axis.title.x = element_blank(),  
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits = c(0, 6)) +
  stat_density(geom = "line", position = "identity", aes(x = pasture1_p50), adjust = 1, kernel = "gaussian")

plot_p75 <- ggplot(mean_ratio_draw, aes(x = pasture1_p75, y = reorder(ecoregion, -pasture1_p75), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.6, scale = 1.2) +  
  labs(title = "75th percentile",
       x = "sensitivity\n(N forest/N pasture)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=6),
        axis.title.x = element_text(size = 10),  
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits = c(0, 6)) +
  stat_density(geom = "line", position = "identity", aes(x = pasture1_p75), kernel = "gaussian")

grid.arrange(plot_p25, plot_p50, plot_p75,  ncol=1)

mean_ratio_draw %>% group_by(ecoregion) %>% summarise(mean=mean(pasture1_mean))

################################################################################
relative_diff <- mean_ratio_draw %>%
  left_join(mean_ratio_col, by = "draw", suffix = c("", "_colombia")) %>% 
  mutate(
    relative_mean_diff = pasture1_mean / pasture1_mean_colombia,
    relative_p25_diff  = pasture1_p25 / pasture1_p25_colombia,
    relative_p50_diff  = pasture1_p50 / pasture1_p50_colombia,
    relative_p75_diff  = pasture1_p75 / pasture1_p75_colombia
  ) %>%
  select(ecoregion, draw, relative_mean_diff, relative_p25_diff, relative_p50_diff, relative_p75_diff)

plot_p25_diff <- ggplot(relative_diff, aes(x = relative_p25_diff, y = reorder(ecoregion, -relative_p25_diff), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.7) +
  labs(title = "25th percentile",
       x = "Relative difference\n(Ecoregion /Colombia)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits = c(0, 7)) +
  stat_density(geom = "line", position = "identity", aes(x = relative_p25_diff), adjust = 1, kernel = "gaussian")

plot_mean_diff <- ggplot(relative_diff, aes(x = relative_mean_diff, y = reorder(ecoregion, -relative_mean_diff), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.6, scale = 1.5) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 0.7) +
  labs(title = "Median",
       x = "Relative difference\n(Ecoregion /Colombia)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=6),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits = c(0, 7)) +
  stat_density(geom = "line", position = "identity", aes(x = relative_mean_diff), adjust = 1, kernel = "gaussian")

plot_p50_diff <- ggplot(relative_diff, aes(x = relative_p50_diff, y = reorder(ecoregion, -relative_p50_diff), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.6, scale = 1.5) +  
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 0.7) +
  labs(title = "Median",
       x = "Relative difference\n(Ecoregion /Colombia)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits = c(0, 7)) +
  stat_density(geom = "line", position = "identity", aes(x = relative_p50_diff), adjust = 1, kernel = "gaussian")

plot_p75_diff <- ggplot(relative_diff, aes(x = relative_p75_diff, y = reorder(ecoregion, -relative_p75_diff), fill = ecoregion)) +
  geom_density_ridges(alpha = 0.6, scale = 1.2) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 0.7) +
  labs(title = "75th percentile",
       x = "Relative difference\n(Ecoregion /Colombia)",
       y = "Ecoregion") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=6),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold")) +
  scale_x_continuous(limits = c(0, 7)) +
  stat_density(geom = "line", position = "identity", aes(x = relative_p75_diff), kernel = "gaussian")

grid.arrange(plot_p25_diff, plot_p50_diff, plot_p75_diff,ncol=1)


grid.arrange(plot_p25, plot_p25_diff, plot_p50, 
             plot_p50_diff, plot_p75, plot_p75_diff,ncol=2)

marta <- ecoregions_predictions[[13]]
marta_forest <- subset(marta, marta$pasture==1)  
marta_pasture <- subset(marta, marta$pasture==0)  

sum(marta_forest$abun__draw_300[marta_forest$scientificName == "Canthidium_sp._01H"], na.rm = TRUE) /
  sum(marta_pasture$abun__draw_300[marta_pasture$scientificName == "Canthidium_sp._01H"], na.rm = TRUE)

################################################################################
# sensitivities by species across Colombia
db_predictions <- readRDS("./species_predictions_2.rds")
# Get the region names from the list names
sp_names <- names(db_predictions)
# calculate the proportional change in abundance for species and region by draw

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
# Apply the function for each species
ratio_draw <- map2(db_predictions, sp_names, abundance_change)
ratio_draw_df <- bind_rows(ratio_draw)
ratio_draw_df$ecoregion <- paste0("Colombia")

mean_ratio <- function(df, df_name) {
  df %>%
    summarise(across(starts_with("ratio_abun__draw_"), 
                     list(mean = ~mean(.x[is.finite(.x)], na.rm = TRUE),
                          p25  = ~quantile(.x[is.finite(.x)], probs = 0.25, na.rm = TRUE),
                          p50  = ~quantile(.x[is.finite(.x)], probs = 0.50, na.rm = TRUE),
                          p75  = ~quantile(.x[is.finite(.x)], probs = 0.75, na.rm = TRUE)))) %>% 
    pivot_longer(cols = starts_with("ratio_abun__draw_"),
                 names_to = c("draw", ".value"),
                 names_pattern = "ratio_abun__draw_(\\d+)_(.*)") %>%
    mutate(ecoregion = df_name)
}
mean_ratio_draw <- mean_ratio(ratio_draw_df, "Colombia")
saveRDS(mean_ratio_draw, "./mean_ratio_col_2.rds")
gc()


