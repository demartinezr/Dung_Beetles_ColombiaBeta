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

