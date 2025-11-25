setwd("C:/Users/PC/Dropbox/CO_DBdata")

library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(raster)
library(terra)
library(patchwork)

# abundance maps for dung beetles in Colombia
#
# abundance predictions for 243 species of dung beetles with 10 iterations in 
# pasture and forest
  db_predictions <- readRDS("C:/Users/PC/Dropbox/CO_DBdata/species_predictions.rds")
# study area ecoregions and grid
  study_area <- readRDS("./SIG/WWF_terrestrial_ecoregions.rds") %>%
    mutate(label_id = 1:n())
  ################################################################################   
  # median and mean abundance for Ateuchus irinus, an endemic species from Magdalena Valley
  irinus <- db_predictions[["Ateuchus_irinus"]]
  irinus <- irinus %>%
    rowwise() %>%
    mutate(median_abundance = median(c_across(starts_with("abun__draw_")), na.rm = TRUE),
           mean_abundance = mean(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
    ungroup()
  
  irinus_pasture <- subset(irinus, irinus$pasture == 1)
  irinus_forest <- subset(irinus, irinus$pasture == 0)
  
  AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  # 3. Reproyectar ambos al sistema Albers
  irinus_forest <- st_transform(irinus_forest, crs = AEAstring)
  irinus_pasture <- st_transform(irinus_pasture, crs = AEAstring)
  study_area <- st_transform(study_area, crs = AEAstring)
  
  colombia <- ne_countries(scale = "medium", returnclass = "sf") %>%
    dplyr::filter(name == "Colombia")
  colombia <- st_transform(colombia, crs = AEAstring)
  dem_col <- rast("F:/Capas/America/dem/elev_raster/raster_elev_AEA.grd")
  dem_df <- as.data.frame(dem_col, xy =TRUE)
  range_irinus <- st_read("F:/Doctorado/Tesis/GBIF/contorno/Ateuchus_irinus.shp")
  
  elevation_dem <- ggplot() +
    geom_sf(data = colombia, fill = "gray", color = "black") +  
    geom_tile(data = dem_df, aes(x = x, y = y, fill = elevation)) + 
    scale_fill_viridis_c() + 
    coord_sf() +
    theme_void() +
    labs(title = "none",
         fill = "none") +
    guides(fill = "none")
  
  #    ggsave("./elevation_dem.jpeg", plot = elevation_dem, width = 2.5, height = 3, units = "in",        
  #           dpi = 300, device = "jpeg")
  
  range_m <- ggplot() +
    geom_sf(data = colombia, fill = "gray90", color = "black") +  
    geom_sf(data = range_irinus, fill = "darkblue", color = "darkblue", alpha = 0.5) +
    coord_sf() +
    theme_void() +
    labs(title = "Distribución de Ateuchus irinus")
  
  #ggsave("./range_m.jpeg", plot = range_m, width = 2.5, height = 3, units = "in",        
  #       dpi = 300, device = "jpeg")          
  
  
  forest_base_irinus <- rast(ext(study_area), 
                             resolution = 2000, 
                             crs = st_crs(study_area)$wkt)
  
  pasture_base_irinus <- rast(ext(study_area), 
                              resolution = 2000, 
                              crs = st_crs(study_area)$wkt)
  
  forest_raster_irinus <- rasterize(vect(irinus_forest), forest_base_irinus, field = "mean_abundance", fun = "mean")
  #    forest_abundance <- global(forest_raster, fun = "sum", na.rm = TRUE)[[1]]
  #    forest_raster <- forest_raster/ forest_abundance
  forest_raster_irinus[is.na(forest_raster_irinus)] <- 0
  forest_raster_irinus[forest_raster_irinus >= 1] <- 1
  #   forest_raster <- log(forest_raster)
  #   forest_min_value <- min(forest_raster[!is.na(forest_raster[]) & forest_raster[] != -Inf], na.rm = TRUE)
  #   forest_raster[forest_raster == -Inf] <- forest_min_value
  forest_raster_masked_irinus <- mask(forest_raster_irinus, vect(study_area))
  forest_irinus_df <- as.data.frame(forest_raster_masked_irinus, xy = TRUE)
  
  p1_irinus <- ggplot(forest_irinus_df) +
    geom_sf(data = colombia, fill = "gray", color = "black") +  
    geom_raster(aes(x = x, y = y, fill = mean)) + 
    scale_fill_viridis_c() + 
    coord_sf() +
    theme_void() +
    guides(fill = "none")+
    theme(legend.position = "none")
  
  pasture_raster_irinus <- rasterize(vect(irinus_pasture), pasture_base_irinus, field = "mean_abundance", fun = "mean")
  pasture_raster_irinus[is.na(pasture_raster_irinus)] <- 0
  #  pasture_raster <- log(pasture_raster)
  #  pasture_min_value <- min(pasture_raster[!is.na(pasture_raster[]) & pasture_raster[] != -Inf], na.rm = TRUE)
  pasture_raster_irinus[pasture_raster_irinus >= 1] <- 1
  #  pasture_raster[pasture_raster == -Inf] <- pasture_min_value
  pasture_raster_masked_irinus <- mask(pasture_raster_irinus, vect(study_area))
  pasture_irinus_df <- as.data.frame(pasture_raster_masked_irinus, xy = TRUE)
  
  # mean abundance map in forest
  p2_irinus <- ggplot(pasture_irinus_df) +
    geom_sf(data = colombia, fill = "gray", color = "black") + 
    geom_raster(aes(x = x, y = y, fill = mean)) + 
    scale_fill_viridis_c() +
    coord_sf() +
    theme_void() +
    theme(legend.position = "none")
  
  fig_1b_a <-p1_irinus + p2_irinus
  ggsave("./fig_1b_a.jpeg", plot = fig_1b_a, width = 8.5, height = 3, units = "in",        
         dpi = 300, device = "jpeg")
  
################################################################################
# median and mean abundance for Ontherus_lunicollis, an endemic species from Magdalena Valley
  lunicollis <- db_predictions[["Ontherus_lunicollis"]]
  lunicollis <- lunicollis %>%
    rowwise() %>%
    mutate(median_abundance = median(c_across(starts_with("abun__draw_")), na.rm = TRUE),
           mean_abundance = mean(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
    ungroup()
  
  lunicollis_pasture <- subset(lunicollis, lunicollis$pasture == 1)
  lunicollis_forest <- subset(lunicollis, lunicollis$pasture == 0)
  
  AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
  
  # 3. Reproyectar ambos al sistema Albers
  lunicollis_forest <- st_transform(lunicollis_forest, crs = AEAstring)
  lunicollis_pasture <- st_transform(lunicollis_pasture, crs = AEAstring)
  study_area <- st_transform(study_area, crs = AEAstring)
  
  colombia <- ne_countries(scale = "medium", returnclass = "sf") %>%
    dplyr::filter(name == "Colombia")
  colombia <- st_transform(colombia, crs = AEAstring)
  dem_col <- rast("F:/Capas/America/dem/elev_raster/raster_elev_AEA.grd")
  dem_df <- as.data.frame(dem_col, xy =TRUE)
  range_lunicollis <- st_read("F:/Doctorado/Tesis/GBIF/contorno/Ontherus_lunicollis.shp")
  
#  elevation_dem <- ggplot() +
#    geom_sf(data = colombia, fill = "gray", color = "black") +  
#    geom_tile(data = dem_df, aes(x = x, y = y, fill = elevation)) + 
#    scale_fill_viridis_c() + 
#    coord_sf() +
#    theme_void() +
#    labs(title = "forest") +
#    guides(fill = "none")
  
#  ggsave("./elevation_dem.jpeg", plot = elevation_dem, width = 2.5, height = 3, units = "in",        
#         dpi = 300, device = "jpeg")
  
  range_m_lunicollis <- ggplot() +
    geom_sf(data = colombia, fill = "gray90", color = "black") +  
    geom_sf(data = range_lunicollis, fill = "darkblue", color = "darkblue", alpha = 0.5) +
    coord_sf() +
    theme_void() +
    labs(title = "Distribución de Ontherus lunicollis")
  
#  ggsave("./range_m.jpeg", plot = range_m, width = 2.5, height = 3, units = "in",        
#         dpi = 300, device = "jpeg")          
  
  
  forest_base_lunicollis <- rast(ext(study_area), 
                      resolution = 2000, 
                      crs = st_crs(study_area)$wkt)
  
  pasture_base_lunicollis <- rast(ext(study_area), 
                       resolution = 2000, 
                       crs = st_crs(study_area)$wkt)
  
  forest_raster_lunicollis <- rasterize(vect(lunicollis_forest), forest_base_lunicollis, field = "mean_abundance", fun = "mean")
  #    forest_abundance <- global(forest_raster, fun = "sum", na.rm = TRUE)[[1]]
  #    forest_raster <- forest_raster/ forest_abundance
  forest_raster_lunicollis[is.na(forest_raster_lunicollis)] <- 0
  forest_raster_lunicollis[forest_raster_lunicollis >= 1] <- 1
  #   forest_raster <- log(forest_raster)
  #   forest_min_value <- min(forest_raster[!is.na(forest_raster[]) & forest_raster[] != -Inf], na.rm = TRUE)
  #   forest_raster[forest_raster == -Inf] <- forest_min_value
  forest_raster_masked_lunicollis <- mask(forest_raster_lunicollis, vect(study_area))
  forest_lunicollis_df <- as.data.frame(forest_raster_masked_lunicollis, xy = TRUE)
  
  p1_lunicollis <- ggplot(forest_lunicollis_df) +
    geom_sf(data = colombia, fill = "gray", color = "black") +  
    geom_raster(aes(x = x, y = y, fill = mean)) + 
    scale_fill_viridis_c() + 
    coord_sf() +
    theme_void() +
    guides(fill = "none")+
    theme(legend.position = "none")
  
  pasture_raster_lunicollis <- rasterize(vect(lunicollis_pasture), pasture_base_lunicollis, field = "mean_abundance", fun = "mean")
  pasture_raster_lunicollis[is.na(pasture_raster_lunicollis)] <- 0
  #  pasture_raster <- log(pasture_raster)
  #  pasture_min_value <- min(pasture_raster[!is.na(pasture_raster[]) & pasture_raster[] != -Inf], na.rm = TRUE)
  pasture_raster_lunicollis[pasture_raster_lunicollis >= 1] <- 1
  #  pasture_raster[pasture_raster == -Inf] <- pasture_min_value
  pasture_raster_masked_lunicollis <- mask(pasture_raster_lunicollis, vect(study_area))
  pasture_lunicollis_df <- as.data.frame(pasture_raster_masked_lunicollis, xy = TRUE)
  
  # mean abundance map in forest
  p2_lunicollis <- ggplot(pasture_lunicollis_df) +
    geom_sf(data = colombia, fill = "gray", color = "black") + 
    geom_raster(aes(x = x, y = y, fill = mean)) + 
    scale_fill_viridis_c() +
    coord_sf() +
    theme_void() +
    theme(legend.position = "none")
  
  fig_1b_b <- p1_lunicollis + p2_lunicollis
  ggsave("./fig_1b_b.jpeg", plot = fig_1b_b, width = 8.5, height = 3, units = "in",        
         dpi = 300, device = "jpeg")
################################################################################   
    # median and mean abundance for Dichotomius_nisus, an endemic species from Magdalena Valley
    nisus <- db_predictions[["Dichotomius_nisus"]]
    nisus <- nisus %>%
      rowwise() %>%
      mutate(median_abundance = median(c_across(starts_with("abun__draw_")), na.rm = TRUE),
             mean_abundance = mean(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
      ungroup()
    
    nisus_pasture <- subset(nisus, nisus$pasture == 1)
    nisus_forest <- subset(nisus, nisus$pasture == 0)
    
    AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    
    # 3. Reproyectar ambos al sistema Albers
    nisus_forest <- st_transform(nisus_forest, crs = AEAstring)
    nisus_pasture <- st_transform(nisus_pasture, crs = AEAstring)
    study_area <- st_transform(study_area, crs = AEAstring)
    
    colombia <- ne_countries(scale = "medium", returnclass = "sf") %>%
      dplyr::filter(name == "Colombia")
    colombia <- st_transform(colombia, crs = AEAstring)
    dem_col <- rast("F:/Capas/America/dem/elev_raster/raster_elev_AEA.grd")
    dem_df <- as.data.frame(dem_col, xy =TRUE)
    range_map <- st_read("F:/Doctorado/Tesis/GBIF/contorno/Dichotomius_nisus.shp")
    
    elevation_dem <- ggplot() +
      geom_sf(data = colombia, fill = "gray", color = "black") +  
      geom_tile(data = dem_df, aes(x = x, y = y, fill = elevation)) + 
      scale_fill_viridis_c() + 
      coord_sf() +
      theme_void() +
      labs(title = "forest",
           fill = "cutoff = 1") +
      guides(fill = "none")
    
#    ggsave("./elevation_dem.jpeg", plot = elevation_dem, width = 2.5, height = 3, units = "in",        
#           dpi = 300, device = "jpeg")
    
    range_m <- ggplot() +
      geom_sf(data = colombia, fill = "gray90", color = "black") +  
      geom_sf(data = range_map, fill = "darkblue", color = "darkblue", alpha = 0.5) +
      coord_sf() +
      theme_void() +
      labs(title = "Distribución de Ontherus lunicollis")
    
#    ggsave("./range_m.jpeg", plot = range_m, width = 2.5, height = 3, units = "in",        
#           dpi = 300, device = "jpeg")          
    
    
    forest_base_nisus <- rast(ext(study_area), 
                        resolution = 2000, 
                        crs = st_crs(study_area)$wkt)
    
    pasture_base_nisus <- rast(ext(study_area), 
                         resolution = 2000, 
                         crs = st_crs(study_area)$wkt)
    
    forest_raster_nisus <- rasterize(vect(nisus_forest), forest_base_nisus, field = "mean_abundance", fun = "mean")
    #    forest_abundance <- global(forest_raster, fun = "sum", na.rm = TRUE)[[1]]
    #    forest_raster <- forest_raster/ forest_abundance
    forest_raster_nisus[is.na(forest_raster_nisus)] <- 0
    forest_raster_nisus[forest_raster_nisus >= 1] <- 1
    #   forest_raster <- log(forest_raster)
    #   forest_min_value <- min(forest_raster[!is.na(forest_raster[]) & forest_raster[] != -Inf], na.rm = TRUE)
    #   forest_raster[forest_raster == -Inf] <- forest_min_value
    forest_raster_masked_nisus <- mask(forest_raster_nisus, vect(study_area))
    forest_nisus_df <- as.data.frame(forest_raster_masked_nisus, xy = TRUE)
    
    p1_nisus <- ggplot(forest_nisus_df) +
      geom_sf(data = colombia, fill = "gray", color = "black") +  
      geom_raster(aes(x = x, y = y, fill = mean)) + 
      scale_fill_viridis_c() + 
      coord_sf() +
      theme_void() +
      guides(fill = "none")+
      theme(legend.position = "none")
    
    pasture_raster_nisus <- rasterize(vect(nisus_pasture), pasture_base_nisus, field = "mean_abundance", fun = "mean")
    pasture_raster_nisus[is.na(pasture_raster_nisus)] <- 0
    #  pasture_raster <- log(pasture_raster)
    #  pasture_min_value <- min(pasture_raster[!is.na(pasture_raster[]) & pasture_raster[] != -Inf], na.rm = TRUE)
    pasture_raster_nisus[pasture_raster_nisus >= 1] <- 1
    #  pasture_raster[pasture_raster == -Inf] <- pasture_min_value
    pasture_raster_masked_nisus <- mask(pasture_raster_nisus, vect(study_area))
    pasture_nisus_df <- as.data.frame(pasture_raster_masked_nisus, xy = TRUE)
    
    # mean abundance map in forest
    p2_nisus <- ggplot(pasture_nisus_df) +
      geom_sf(data = colombia, fill = "gray", color = "black") + 
      geom_raster(aes(x = x, y = y, fill = mean)) + 
      scale_fill_viridis_c() +
      coord_sf() +
      theme_void() +
      theme(legend.position = "none")
    
    fig_1b_c <- p1_nisus + p2_nisus
    ggsave("./fig_1b_c.jpeg", plot = fig_1b_c, width = 8.5, height = 3, units = "in",        
           dpi = 300, device = "jpeg")
    
    
    fig_1_pred <- fig_1b_a / fig_1b_c
    
    ggsave("./fig_1_predA.jpeg", plot = fig_1_pred, width = 4.5, height = 6, units = "in",        
           dpi = 300, device = "jpeg")
################################################################################      
################################################################################
###### cutoff_3 or 10 change filter raster
    
    forest_raster <- rasterize(vect(cracicus_forest), forest_base, field = "mean_abundance", fun = "mean")
    forest_raster[is.na(forest_raster)] <- 0
    forest_raster[forest_raster >= 10] <- 1
    forest_raster_masked <- mask(forest_raster, vect(study_area))
    forest_df <- as.data.frame(forest_raster_masked, xy = TRUE)

    p1 <- ggplot(forest_df) +
      geom_sf(data = colombia, fill = "gray", color = "black") +  
      geom_raster(aes(x = x, y = y, fill = mean)) + 
      scale_fill_viridis_c() + 
      coord_sf() +
      theme_void() +
      labs(title = "forest",
           fill = "cutoff = 10") +
      guides(fill = "none")
    
    pasture_raster <- rasterize(vect(cracicus_pasture), pasture_base, field = "mean_abundance", fun = "mean")
    pasture_raster[is.na(pasture_raster)] <- 0
    pasture_raster[pasture_raster >= 10] <- 1
    pasture_raster_masked <- mask(pasture_raster, vect(study_area))
    pasture_df <- as.data.frame(pasture_raster_masked, xy = TRUE)
    
    # mean abundance map in forest
    p2 <- ggplot(pasture_df) +
      geom_sf(data = colombia, fill = "gray", color = "black") + 
      geom_raster(aes(x = x, y = y, fill = mean)) + 
      scale_fill_viridis_c() +
      coord_sf() +
      theme_void() +
      labs(title = "Pasture",
           fill = "cutoff = 10")
    
    fig_1b_3_cut10 <- grid.arrange(p1, p2, ncol=2, widths =c(0.9,1))
    saveRDS(fig_1b_3_cut10, "./pictures/fig_1b_3_cut10.rds")
    
    
    ggsave("./fig_1b_3_cut10.jpeg", plot = fig_1b_3_cut10, width = 8.5, height = 3, units = "in",        
           dpi = 300, device = "jpeg")

  # points to grid 2km for A. cracicus in pasture 
  
    cracicus_grid_pasture <- cracicus_pasture  %>%
      mutate(geometry = st_buffer(geometry, dist = grid_size / sqrt(2), endCapStyle = "SQUARE")) %>%
      st_as_sf()
  
    p3 <- ggplot() +
      geom_sf(data = cracicus_grid_pasture, aes(fill = log(median_abundance)), color = NA) +  
      scale_fill_viridis_c(option = "C") +
      theme_classic() +
      labs(title = "Ateuchus cracicus in pasture",
           fill = "Median abundance")
    
    p4 <- ggplot() +
      geom_sf(data = cracicus_grid_pasture, aes(fill = log(mean_abundance)), color = NA) +  
      scale_fill_viridis_c(option = "C", limits = c(-2, 2)) +
      theme_classic() +
      labs(title = "Ateuchus cracicus in pasture",
           fill = "Mean abundance")
    
    grid.arrange(p1, p2, p3, p4, ncol=2) 
    grid.arrange(p2, p4, ncol=2) 
#
# Using a grid 2km and join predictions points for abundance maps
  # join abundance to grid in forest
  cracicus_col_forest <- grid_2km %>%
    st_join(cracicus_forest, left = TRUE)
  #cracicus_col_forest <- cracicus_col_forest %>%
    #mutate(mean_abundance = ifelse(is.na(mean_abundance), 0, mean_abundance))
  # abundance map in forest
  p5 <- ggplot(cracicus_col_forest) +
          geom_sf(aes(fill = mean_abundance), color = NA) + 
          scale_fill_viridis_c(option = "C") +
          theme_classic() +
          labs(title = "Ateuchus cracicus forest",,
               fill = "Mean abundance")
  # join abundance to grid in pasture      
  cracicus_col_pasture <- grid_2km %>%
    st_join(cracicus_pasture, left = TRUE)
    #cracicus_col_pasture <- cracicus_col_pasture %>%
      #mutate(mean_abundance = ifelse(is.na(mean_abundance), 0, mean_abundance))
  # abundance map in pasture
  p6 <- ggplot(cracicus_col_pasture) +
          geom_sf(aes(fill = mean_abundance), color = NA) + 
          scale_fill_viridis_c(option = "C") +
          theme_classic() +
          labs(title = "Ateuchus cracicus in pasture",
               fill = "Mean abundance")
        
  grid.arrange(p5, p6, ncol=2)
#
# median and mean abundance for Digitonthophagus gazella an introduced species
#
  gazella <- db_predictions[["Digitonthophagus_gazella"]]
  gazella <- gazella %>%
    rowwise() %>%
    mutate(median_abundance = median(c_across(starts_with("abun__draw_")), na.rm = TRUE),
           mean_abundance = mean(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
    ungroup()
  
    gazella_pasture <- subset(gazella, gazella$pasture == 0)
    gazella_forest <- subset(gazella, gazella$pasture == 1)
  
  # points to grid 2km for D. gazella in forest
    grid_size <- 2000
    gazella_grid_forest <- gazella_forest %>%
      mutate(geometry = st_buffer(geometry, dist = grid_size / sqrt(2), endCapStyle = "SQUARE")) %>%
      st_as_sf()
    # median abundance map in forest
    p7 <- ggplot() +
            geom_sf(data = gazella_grid_forest, aes(fill = median_abundance), color = "NA") +  
            scale_fill_viridis_c(option = "C") +
            theme_classic() +
            labs(title = "Digitonthophagus gazella in forest",
                 fill = "Median abundance")
    # mean abundance map in forest
    p8 <- ggplot() +
      geom_sf(data = gazella_grid_forest, aes(fill = mean_abundance), color = "NA") +  
      scale_fill_viridis_c(option = "C") +
      theme_classic() +
      labs(title = "Digitonthohpagus gazella in forest",
           fill = "Mean abundance")
    grid.arrange(p7, p8, ncol=2)
    
  # points to grid 2km for D. gazella in pasture 
    gazella_grid_pasture <- gazella_pasture  %>%
      mutate(geometry = st_buffer(geometry, dist = grid_size / sqrt(2), endCapStyle = "SQUARE")) %>%
      st_as_sf()
    # median abundance map in pasture
    p9 <- ggplot() +
      geom_sf(data = gazella_grid_pasture, aes(fill = median_abundance), color = NA) +  
      scale_fill_viridis_c(option = "C") +
      theme_classic() +
      labs(title = "Digitonthophagus gazella in pasture",
           fill = "Median abundance")
    # mean abundance map in pasture
    p10 <- ggplot() +
      geom_sf(data = gazella_grid_pasture, aes(fill = mean_abundance), color = NA) +  
      scale_fill_viridis_c(option = "C") +
      theme_classic() +
      labs(title = "Digitonthophagus gazella in pasture",
           fill = "Mean abundance")
  
    grid.arrange(p7, p8, p9, p10, ncol=2)  
  #
  # join abundance to grid in forest
  gazella_col_forest <- grid_2km %>%
    st_join(gazella_forest, left = TRUE)
  #gazella_col_forest <- gazella_col_forest %>%
    #mutate(median_abundance = ifelse(is.na(median_abundance), 0, median_abundance))
  # abundance map in forest
  p11 <- ggplot(gazella_col_forest) +
          geom_sf(aes(fill = mean_abundance), color = NA) + 
          scale_fill_viridis_c(option = "C")  +
          theme_classic() +
          labs(title = "Digitonthophagus gazella forest")
  # join abundance to grid in pasture     
  gazella_col_pasture <- grid_2km %>%
    st_join(gazella_pasture, left = TRUE)
  #gazella_col_pasture <- gazella_col_pasture %>%
    #mutate(median_abundance = ifelse(is.na(median_abundance), 0, median_abundance))
  # abundance map in pasture  
   p12 <- ggplot(gazella_col_pasture) +
           geom_sf(aes(fill = mean_abundance), color = NA) + 
           scale_fill_viridis_c(option = "C") +
           theme_classic() +
           labs(title = "Digitonthophagus gazella pasture")
#
# median and mean abundancia for Dichotomius andresi, an endemic species of Magdalena Valley
#
  andresi <- db_predictions[["Dichotomius_andresi"]]
  andresi <- andresi %>%
      rowwise() %>%
      mutate(median_abundance = median(c_across(starts_with("abun__draw_")), na.rm = TRUE),
             mean_abundance = mean(c_across(starts_with("abun__draw_")), na.rm = TRUE)) %>%
      ungroup()
    
  andresi_pasture <- subset(andresi, andresi$pasture == 0)
  andresi_forest <- subset(andresi, andresi$pasture == 1)
    
# points to grid 2km for D. andresi in forest
    grid_size <- 2000
    andresi_grid_forest <- andresi_forest %>%
      mutate(geometry = st_buffer(geometry, dist = grid_size / sqrt(2), endCapStyle = "SQUARE")) %>%
      st_as_sf()
    # median abundance map in forest
    p13 <- ggplot() +
      geom_sf(data = andresi_grid_forest, aes(fill = log10(median_abundance)), color = "NA") +  
      scale_fill_viridis_c(option = "C") +
      theme_classic() +
      labs(title = "Dichotomius andresi in forest",
           fill = "Median abundance")
    # mean abundance map in forest
    p14 <- ggplot() +
      geom_sf(data = andresi_grid_forest, aes(fill = log10(mean_abundance)), color = "NA") +  
      scale_fill_viridis_c(option = "C") +
      theme_classic() +
      labs(title = "Dichotomius andresi in forest",
           fill = "Mean abundance")
    grid.arrange(p13, p14, ncol=2)
    
    # points to grid 2km for D. andresi in pasture 
    andresi_grid_pasture <- andresi_pasture  %>%
      mutate(geometry = st_buffer(geometry, dist = grid_size / sqrt(2), endCapStyle = "SQUARE")) %>%
      st_as_sf()
    # median abundance map in pasture
    p15 <- ggplot() +
      geom_sf(data = andresi_grid_pasture, aes(fill = log10(median_abundance)), color = NA) +  
      scale_fill_viridis_c(option = "C") +
      theme_classic() +
      labs(title = "Dichotomius andresi in pasture",
           fill = "Median abundance")
    # mean abundance map in pasture
    p16 <- ggplot() +
      geom_sf(data = andresi_grid_pasture, aes(fill = log10(mean_abundance)), color = NA) +  
      scale_fill_viridis_c(option = "C") +
      theme_classic() +
      labs(title = "Dichotomius andresi in pasture",
           fill = "Mean abundance")
    
    grid.arrange(p13, p14, p15, p16, ncol=2)  
    #
  # join abundance to grid in forest
    andresi_col_forest <- grid_2km %>%
      st_join(andresi_forest, left = TRUE)
    andresi_col_forest <- andresi_col_forest %>%
    mutate(mean_abundance = ifelse(is.na(mean_abundance), -2, mean_abundance))
    # abundance map in forest
    p17 <- ggplot(andresi_col_forest) +
      geom_sf(aes(fill = log10(mean_abundance)), color = NA) + 
      scale_fill_viridis_c(option = "C")  +
      theme_classic() +
      labs(title = "Dichotomius andresi forest")
  # join abundance to grid in pasture     
    andresi_col_pasture <- grid_2km %>%
      st_join(andresi_pasture, left = TRUE)
    andresi_col_pasture <- andresi_col_pasture %>%
    mutate(mean_abundance = ifelse(is.na(mean_abundance), -2, mean_abundance))
    # abundance map in pasture  
    p18 <- ggplot(andresi_col_pasture) +
      geom_sf(aes(fill = log10(mean_abundance)), color = NA) + 
      scale_fill_viridis_c(option = "C") +
      theme_classic() +
      labs(title = "Dichotomius andresi pasture")
    grid.arrange(p17, p18, ncol=2)  
    
