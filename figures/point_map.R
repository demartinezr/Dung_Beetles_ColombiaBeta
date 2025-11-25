setwd("C:/Users/PC/Dropbox/CO_DBdata")

library(sf)
library(terra)
library(ggplot2)
library(tidyverse)
library(ggspatial)
library(ggnewscale)
library(RColorBrewer)
library(cowplot) 


hillshade <- rast("F:/Capas/Hillshade/World_e-Atlas-UCSD_SRTM30-plus_v8_Hillshading.tiff")
countries <- st_read("F:/Capas/America/countries/America_countriesC.shp")
ocean <- st_read("F:/Capas/World/ne_10m_ocean/ne_10m_ocean.shp")
points <- readRDS("C:/Users/PC/Dropbox/CO_DBdata/survey/all_pts.RDS")
points <- points[!is.na(points$beetles), ]
points <- st_as_sf(points, coords = c("lon", "lat"), crs = 4326)
points <- points[points$habitat !="Sy",]
points <- points[points$habitat !="PALM",]

countries_crop <- st_crop(countries, ext(-79, -64, -4.6, 13.1))

hillshade_crop <- crop(hillshade, ext(-79, -64, -4.6, 13.1))
hillshade_df <- as.data.frame(hillshade_crop, xy = TRUE, na.rm = TRUE)
names(hillshade_df)[3] <- "elevation"

WWF_ecoregions <- readRDS("./SIG/WWF_terrestrial_ecoregions.rds") %>%
  mutate(label_id = 1:n())

# Fix typos
WWF_ecoregions$ECO_NAME[11] <- "Caquetá moist forests"

#Region full names
regions <- c("Eastern Cordillera real montane forests",
             "Magdalena-Urabá moist forests",
             "Apure-Villavicencio dry forests",
             "Napo moist forests",
             "Cordillera Oriental montane forests",
             "Llanos",
             "Santa Marta montane forests",
             "Magdalena Valley montane forests",
             "Northern Andean páramo",
             "Caquetá moist forests",
             "Magdalena Valley dry forests",
             "Cauca Valley montane forests",
             "Northwestern Andean montane forests",
             "Colombia")

WWF_ecoregions$ECO_NAME <- factor(WWF_ecoregions$ECO_NAME, levels = regions)

# Region short names
reg_short <- c("CC montane",
               "Magdalena-Urabá moist",
               "Villavicencio dry",
               "Napo moist",
               "EC montane",
               "Llanos",
               "Santa Marta montane",
               "Magdalena montane",
               "Andean páramo",
               "Caquetá moist",
               "Magdalena dry",
               "Cauca montane",
               "WC montane",
               "Colombia")
label_vec <- setNames(reg_short, regions)

# Create a color palette function based on RdBu (red-blue)
colors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(13))
col_vec <- setNames(colors, regions[-14])

map_1 <- ggplot() +
  # Hillshade
  geom_raster(data = hillshade_df, aes(x = x, y = y, fill = elevation), show.legend = FALSE) +
  scale_fill_gradient(low = "gray20", high = "gray95", guide = "none") +
  # Océanos
  geom_sf(data = ocean, fill = "lightblue", color = NA,  alpha = 0.7) +
  # Ecoregiones
  new_scale_fill() +
  geom_sf(data = WWF_ecoregions, aes(fill = ECO_NAME), color = "black", size = 0.3, alpha = 0.6) +
  #scale_fill_manual(values = col_vec) +  
  scale_fill_manual(values = col_vec, labels = label_vec) +
  # Limites de país
  geom_sf(data = countries_crop, fill = NA, color = "black", size = 0.6) +
  geom_sf(data = points, fill = "white", color = "black", size = 3, shape = 21) +
  annotation_scale(location = "bl", width_hint = 0.4, line_width = 1, text_cex = 0.8) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         style = north_arrow_fancy_orienteering(fill = c("black", "white")),
                         height = unit(1.3, "cm"), width = unit(1.3, "cm")) +
  coord_sf(xlim = c(-79, -64), ylim = c(-4.6, 13.1), expand = FALSE) +
  annotate("text",
           x = c(-79, -72, -64),
           y = rep(13.2, 3),
           label = c("80°W", "72°W", "66°W"),
           size = 2.5, hjust = 0.5) +
  theme_bw() +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.99, 0.99),         # x (izq) y (arriba)
    legend.justification = c(1, 1),                 # esquina superior izquierda
    legend.background = element_rect(fill = "gray", color = "black", linewidth = 0.5),
    legend.key.size = unit(1, "lines"),
    legend.title = element_blank(),
    legend.text = element_text(size = 11, hjust = 0, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )  +
  labs(color = "Ecoregions")

ggsave("mapa_ecoregiones10.jpeg", map_1, width = 5.5, height = 6.5, units = "in", dpi = 300)



###############################################################################
# sensitivity across posterior
ecoregions_all <- as.data.frame(readRDS("./ecoregions_100.rds"))

library(purrr)
library(dplyr)
library(tidyr)
library(stringr)
library(magick)

ecoregions_all <- as.data.frame(readRDS("./ecoregions_100.rds"))
#ecoregion_sp <- split(ecoregions_all, ecoregions_all$scientificName)
#saveRDS(ecoregion_sp, "./diversity_loss/sp_ecoregion.rds")
ecoregion_sp <- readRDS("./diversity_loss/sp_ecoregion.rds")
# Get the region names from the list names
sp_names <- names(ecoregion_sp)
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

###############################################################################
library(magick)
library(ggplot2)
library(cowplot)
fig_a <- image_read("./fig_1A.jpeg")
fig_b <- image_read("./fig_1C_A.jpg")
fig_c <- image_read("./fig_1C.jpg")
fig_d <- image_read("./fig_1d_rot.jpg")
fig_A <- ggdraw() + draw_image(fig_a)
fig_B <- ggdraw() + draw_image(fig_b)
fig_C <- ggdraw() + draw_image(fig_c)
fig_D <- ggdraw() + draw_image(fig_d)

final_plot_cowplot <- plot_grid(
  plot_grid(fig_A, fig_B, ncol = 2, labels = c("a", "b"), label_size = 12, label_fontface = "bold"),
  ncol = 1)
final_plot_cowplotA <- plot_grid(
  plot_grid(fig_A, fig_C, ncol = 2, labels = c("a", "b"), label_size = 12, label_fontface = "bold"),
  ncol = 1)


ggsave("fig_1.jpg", final_plot_cowplotA, width = 8.5, height = 5.4, units = "in", dpi = 300)

# Crear figura vacía de fondo (tamaño carta)
final_plot <- ggdraw() +
  draw_plot(fig_A, x = 0, y = 0.5, width = 2/3, height = 0.5) +
  draw_label("A", x = 0.02, y = 0.97, fontface = "bold", size = 14) +
  draw_plot(fig_B, x = 2/3, y = 0.5, width = 1/3, height = 0.5) +
  draw_label("B", x = 0.68, y = 0.97, fontface = "bold", size = 14) +
  draw_plot(fig_C, x = 0, y = 0, width = 0.5, height = 0.5) +
  draw_label("C", x = 0.02, y = 0.47, fontface = "bold", size = 14) +
  draw_plot(fig_D, x = 0.5, y = 0, width = 0.5, height = 0.5) +
  draw_label("D", x = 0.52, y = 0.47, fontface = "bold", size = 14)
  
ggsave("final_plot.jpg", final_plot, width = 8.5, height = 11, units = "in", dpi = 300)
