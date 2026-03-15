setwd("C:/Users/PC/Dropbox/CO_DBdata")

# sensitivity across posterior
library(sf)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(purrr)
library(stringr)
library(magick)

ecoregions_all <- as.data.frame(readRDS("./ecoregions_100.rds"))
#ecoregion_sp <- split(ecoregions_all, ecoregions_all$scientificName)
#saveRDS(ecoregion_sp, "./diversity_loss/sp_ecoregion.rds")
ecoregion_sp <- readRDS("./diversity_loss/sp_ecoregion.rds")
# Get the region names from the list names
sp_names <- names(ecoregion_sp)

df <- as.data.frame(ecoregion_sp$Anisocanthon_villosus)
df <- st_drop_geometry(df)

df_long <- df %>%
  pivot_longer(
    cols = starts_with("abun__draw_"),
    names_to = "draw",
    values_to = "abundance"
  ) %>%
  mutate(
    draw = as.integer(gsub("abun__draw_", "", draw))
  )

df_ratio <- df_long %>%
  select(scientificName, lon, lat, ecoregions, draw, pasture, abundance) %>%
  pivot_wider(
    names_from = pasture,
    values_from = abundance,
    names_prefix = "pasture_"
  ) %>%
  mutate(
    ratio_fp = pasture_0 / pasture_1
  )

ratio_2km_summary <- df_ratio %>%
  filter(is.finite(ratio_fp)) %>%
  group_by(scientificName, draw) %>%
  summarise(
    ratio_mean   = mean(ratio_fp),
    ratio_median = median(ratio_fp),
    n_cells      = n(),
    .groups = "drop"
  )

############################## 2 x 2 km

ratio_2km_sp <- function(df) {
  
  df <- df %>%
    st_drop_geometry() %>%
    as.data.frame()
  
  df_long <- df %>%
    pivot_longer(
      cols = starts_with("abun__draw_"),
      names_to = "draw",
      values_to = "abundance"
    ) %>%
    mutate(
      draw = as.integer(gsub("abun__draw_", "", draw))
    )
  
  df_ratio <- df_long %>%
    select(scientificName, lon, lat, ecoregions, draw, pasture, abundance) %>%
    pivot_wider(
      names_from = pasture,
      values_from = abundance,
      names_prefix = "pasture_"
    ) %>%
    mutate(
      ratio_fp = pasture_0 / pasture_1,
      ratio_fp = ifelse(
        is.infinite(ratio_fp) | is.nan(ratio_fp),
        NA_real_,
        ratio_fp
      )
    )
  
  ratio_2km_summary <- df_ratio %>%
    group_by(scientificName, draw) %>%
    summarise(
      ratio_mean   = mean(ratio_fp, na.rm = TRUE),
      ratio_median = median(ratio_fp, na.rm = TRUE),
      n_cells      = sum(!is.na(ratio_fp)),
      .groups = "drop"
    )
  
  return(ratio_2km_summary)
}

ratio_2km_list <- map(
  ecoregion_sp,
  ratio_2km_sp
)

#saveRDS(ratio_2km_list, "./Analysis/ratio_2km_list.rds")

ratio_2km_list <- readRDS("./Analysis/ratio_2km_list.rds")
ratio_2km_all <- bind_rows(ratio_2km_list)

median_2km <- ratio_2km_all %>%
  group_by(scientificName) %>%
  summarise(
    median_ratio = mean(ratio_median, na.rm = TRUE),
    .groups = "drop"
  )

fig_2km <- ggplot(median_2km, aes(x = log(median_ratio))) +
  geom_histogram(binwidth = 0.2, fill = "grey90", color = "black", alpha = 1) +
  geom_vline(xintercept = 0, linetype = "solid", color = "red", linewidth = 1) +
  #  geom_vline(xintercept = c(-2.5, -0.1, 6.4), 
  #             linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_y_continuous(limits=c(0, 15)) +
  scale_x_continuous(limits= c(-5, 12), 
                     breaks = log(c(0.1, 1, 10, 1000)),
                     labels = c("10x\npasture", "0", "10x", "1000x\nforest")) +
  labs(title = "Local (2 x 2 km)", y = "Frecuency", x = " ") +
  theme_bw() +
  theme(plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 11)) 

# ----------- ecoregions
ratio_eco_sp <- function(df) {
  
  df <- df %>%
    st_drop_geometry() %>%
    as.data.frame()
  
  df_long <- df %>%
    pivot_longer(
      cols = starts_with("abun__draw_"),
      names_to = "draw",
      values_to = "abundance"
    ) %>%
    mutate(
      draw = as.integer(gsub("abun__draw_", "", draw))
    )
  
  df_ratio <- df_long %>%
    group_by(scientificName, ecoregions, draw, pasture) %>%
    summarise(
      abundance = sum(abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = pasture,
      values_from = abundance,
      names_prefix = "pasture_"
    ) %>%
    mutate(
      ratio_fp = pasture_0 / pasture_1,
      ratio_fp = ifelse(
        is.infinite(ratio_fp) | is.nan(ratio_fp),
        NA_real_,
        ratio_fp
      )
    )
  
  ratio_eco_summary <- df_ratio %>%
    group_by(scientificName, draw) %>%
    summarise(
      ratio_mean   = mean(ratio_fp, na.rm = TRUE),
      ratio_median = median(ratio_fp, na.rm = TRUE),
      n_cells      = sum(!is.na(ratio_fp)),
      .groups = "drop"
    )
  
  return(ratio_eco_summary)
}

ratio_eco_list <- map(
  ecoregion_sp,
  ratio_eco_sp
)

#saveRDS(ratio_eco_list, "./Analysis/ratio_eco_list.rds")

ratio_eco_list <- readRDS("./Analysis/ratio_eco_list.rds")
ratio_eco_all <- bind_rows(ratio_eco_list)

median_eco <- ratio_eco_all %>%
  group_by(scientificName) %>%
  summarise(
    median_ratio = mean(ratio_median, na.rm = TRUE),
    .groups = "drop"
  )

fig_eco <- ggplot(median_eco, aes(x = log(median_ratio))) +
  geom_histogram(binwidth = 0.2, fill = "grey90", color = "black", alpha = 1) +
  geom_vline(xintercept = 0, linetype = "solid", color = "red", linewidth = 1) +
  #  geom_vline(xintercept = c(-2.5, -0.1, 6.4), 
  #             linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_y_continuous(limits=c(0, 15)) +
  scale_x_continuous(limits= c(-5, 15), 
                     breaks = log(c(0.1, 1, 10, 1000)),
                     labels = c("10×\npasture", "0", "10×", "1000×\nforest")) +
  labs(title="Ecoregion", y = NULL, x= "Sensitivity") +
  theme_bw() +
  theme(plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 11)) 

#--------------- Near national
ratio_near_sp <- function(df) {
  
  df <- df %>%
    st_drop_geometry() %>%
    as.data.frame()
  
  df_long <- df %>%
    pivot_longer(
      cols = starts_with("abun__draw_"),
      names_to = "draw",
      values_to = "abundance"
    ) %>%
    mutate(
      draw = as.integer(gsub("abun__draw_", "", draw))
    )
  
  df_ratio <- df_long %>%
    group_by(scientificName, draw, pasture) %>%
    summarise(
      abundance = sum(abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = pasture,
      values_from = abundance,
      names_prefix = "pasture_"
    ) %>%
    mutate(
      ratio_fp = pasture_0 / pasture_1,
      ratio_fp = ifelse(
        is.infinite(ratio_fp) | is.nan(ratio_fp),
        NA_real_,
        ratio_fp
      )
    )
  
  ratio_near_summary <- df_ratio %>%
    group_by(scientificName, draw) %>%
    summarise(
      ratio_mean   = mean(ratio_fp, na.rm = TRUE),
      ratio_median = median(ratio_fp, na.rm = TRUE),
      n_cells      = sum(!is.na(ratio_fp)),
      .groups = "drop"
    )
  
  return(ratio_near_summary)
}

ratio_near_list <- map(
  ecoregion_sp,
  ratio_near_sp
)

#saveRDS(ratio_near_list, "./Analysis/ratio_near_list.rds")

ratio_near_list <- readRDS("./Analysis/ratio_near_list.rds")
ratio_near_all <- bind_rows(ratio_near_list)

median_near <- ratio_near_all %>%
  group_by(scientificName) %>%
  summarise(
    median_ratio = mean(ratio_median, na.rm = TRUE),
    .groups = "drop"
  )

fig_near <- ggplot(median_near, aes(x = log(median_ratio))) +
  geom_histogram(binwidth = 0.2, fill = "grey90", color = "black", alpha = 1) +
  geom_vline(xintercept = 0, linetype = "solid", color = "red", linewidth = 1) +
  #  geom_vline(xintercept = c(-2.5, -0.1, 6.4), 
  #             linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_y_continuous(limits=c(0, 15)) +
  scale_x_continuous(limits= c(-5, 15),
                     breaks = log(c(0.1, 1, 10, 1000)),
                     labels = c("10x\npasture", "0", "10x", "1000x\nforest")) +
  labs(title = "Near national", y = NULL, x = " ") +
  theme_bw() +
  theme(plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 11)) 

#---------------
library(patchwork)
fig_1c <- (fig_2km + fig_eco + fig_near)

ggsave("./fig_1c.jpeg", plot = fig_1c, width = 8.5, height = 3, units = "in",        
       dpi = 300, device = "jpeg")

##################################

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
ratio_draw <- map2(ecoregion_sp, sp_names, abundance_change)
#saveRDS(ratio_draw, "./diversity_loss/sp_ratio.rds")

ratio_draw <- readRDS("./diversity_loss/sp_ratio.rds")
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
fig_1d <- ggplot(mean_ratio_draw, aes(x = log(median_ratio))) +
  geom_histogram(binwidth = 0.2, fill = "grey90", color = "black", alpha = 1) +
  geom_vline(xintercept = 0, linetype = "solid", color = "red", linewidth = 1) +
  #  geom_vline(xintercept = c(-2.5, -0.1, 6.4), 
  #             linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_x_continuous(
    name = "Sensitivity",
    breaks = log(c(0.1, 1, 1000)),
    labels = c("10x\npasture", "0", "1000x\nforest")) +
  labs(y = "Frequency") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 11)) 

ggsave("./fig_1d.jpeg", plot = fig_1d, width = 5, height = 2.5, units = "in",        
       dpi = 300, device = "jpeg")

img <- image_read("fig_1d.jpeg")
img_rotated <- image_rotate(img, -90)
image_browse(img_rotated)
image_write(img_rotated, path = "fig_1d_rot.jpg", format = "jpeg")

especies_pasture_dominant <- mean_ratio_draw %>%
  filter(median_ratio < 1)


ratio_col <- readRDS("./decline/ratio_col.rds")
ratio_col$lat_lon <- paste(ratio_col$lat, "_", ratio_col$lon)


paramo <- subset(ratio_col_df, ratio_col$ecoregions =="Northern Andean páramo")

paramo_long <- paramo %>%
  pivot_longer(
    cols = starts_with("ratio__draw_"),
    names_to = "draw",
    names_prefix = "ratio__draw_",
    values_to = "ratio"
  ) %>%
  mutate(draw = as.integer(draw))


sensitivity_draw <- paramo_long %>%
  filter(
    is.finite(ratio)  # elimina NA, NaN, Inf, -Inf
  ) %>%
  group_by(draw, scientificName) %>%
  summarise(
    mean = mean(ratio, na.rm = TRUE),
    median = median(ratio, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

sensitivity_pointwise <- paramo_long %>%
  filter(
    is.finite(ratio)  # elimina NA, NaN, Inf, -Inf
  ) %>%
  group_by(scientificName) %>%
  summarise(
    mean = mean(ratio, na.rm = TRUE),
    median = median(ratio, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

paramo_pasture_pointwise <- sensitivity_pointwise %>%
  filter(median < 1)