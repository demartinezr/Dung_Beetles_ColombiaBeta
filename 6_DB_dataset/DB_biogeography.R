setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata/SIG/inputs")

library(dplyr)
library(ggplot2)
library(rgdal)
library(sf)
library(readxl)


# Dung beetles data set project all projects Western cordillera, llanos, Paramos, beta diversity
  db <- read_excel("C:/Users/Dell-PC/Dropbox/CO_dbdata/abundance/Scarabaeinae_database_2024.xlsx", sheet = "Scarabaeinae_database_2024")
  db$scientificName1 <- gsub("_", " ", db$scientificName)
  db$decimalLongitude <- as.numeric(db$lon_all_points)
  db$decimalLatitude <- as.numeric(db$lat_all_points)
  db$coordinates <- paste(db$lon_all_points, db$lat_all_points, sep = "_")
#
# All sites visited in the projects: Western cordillera, llanos, Paramos, beta diversity
  points <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/all_pts.rds")
  pts_add <- data.frame(point = "TAP1", lat = 5.997912968, lon = -73.220879,
            site = "El Taladro", cluster = "cluster_TAP_1", birds = 0, beetles = 1,
    habitat = "Pasture", natural = 0, paramo = 0, pasture = 1,other = 0,mixed_cluster = 0,
    elev_ALOS = 2350, posix1 = NA, posix2 = NA, posix3 = NA, posix4 = NA, obs1 = NA,
    obs2 = NA, obs3 = NA, obs4 = NA, oday1 = NA, oday2 = NA, oday3 = NA, oday4 = NA,
    hps1 = NA, hps2 = NA, hps3 = NA, hps4 = NA, nv = 4, subregion = "subregion_bogotacito",
    stringsAsFactors = FALSE
  )
  points <- rbind(points, pts_add)
  #
# unifying traps codes between datasets
  db$point[db$point %in% c("SEP_1", "SEP_2", "SEP_3")] <- 
  gsub("_", "", db$point[db$point %in% c("SEP_1", "SEP_2", "SEP_3")])
  #
# obtain site and cluster codes from all pts to dung beetles dataset
  cluster_codes <- points[, c("point", "site", "cluster")]
  db <- merge(db, cluster_codes, by="point", all.x=TRUE)
  #sort(unique(db$point))[!(sort(unique(db$point)) %in% 
  #                                 sort(unique(points$point)))]
  #
# dung beetle biogeography format by Jacob Socolar, adapted code by DEMR with new data
dist <- read_excel("C:/Users/Dell-PC/Dropbox/CO_DBdata/DB_distri-20-07-21.xlsx", sheet="DB_distri-20-07-21")

`%ni%` <- Negate(`%in%`)

sum(db$point %ni% points$point)
unique(db$point[db$point %ni% points$point])

sum(points$beetles == 1) # total traps visited for dung beetles and birds
# unique(points$point[points$point %ni% db$point]) # traps test codes

source("D:/Doctorado/Tesis/Topographic_Col/unidades_topograficas.R")

points_sf <- st_as_sf(points, coords = c("lon", "lat"), crs = 4326)
points_sf$dist_slope <- points_sf$dist_region <- NA

points_sf$dist_region[st_intersects(points_sf, snsm, sparse = F)] <- "SNS Marta"
points_sf$dist_region[!st_intersects(points_sf, snsm, sparse = F) &
                        ((!st_intersects(points_sf, amazon_orinoco, sparse = F)) |
                           (points_sf$elev_ALOS > 500))] <- "Andean"
points_sf$dist_region[grepl("leguizamo", points_sf$subregion)] <- "Amazon"
points_sf$dist_region[grepl("chiribiquete", points_sf$subregion)] <- "Amazon"
points_sf$dist_region[grepl("guaviare", points_sf$subregion)] <- "Amazon"
points_sf$dist_region[grepl("llanos", points_sf$subregion)] <- "Llanos"



points_sf$dist_slope[st_intersects(points_sf, amazon_orinoco, sparse = F) |
                       st_intersects(points_sf, catatumbo, sparse = F)] <- "EC: Eastern"
points_sf$dist_slope[st_intersects(points_sf, magdalena_east, sparse = F)] <- "EC: Western"
points_sf$dist_slope[st_intersects(points_sf, magdalena_west, sparse = F)] <- "CC: Eastern"
points_sf$dist_slope[st_intersects(points_sf, cauca_east, sparse = F)] <- "CC: Western"
points_sf$dist_slope[st_intersects(points_sf, cauca_west, sparse = F)] <- "WC: Eastern"
points_sf$dist_slope[st_intersects(points_sf, pacific, sparse = F)] <- "WC: Western"
points_sf$dist_slope[grep("PUP", points_sf$point)] <- "CC: Western slope"

db1 <- db[!is.na(db$abundance), ]
db1 <- db1[db$abundance > 0, ]
db1 <- db1[!(db1$scientificName=="NA"),]

counter1 <- 0
counter2 <- 0
region_mismatch <- slope_mismatch <- vector()

for(i in 1:nrow(db1)) {
  point <- db1$point[i]
  species <- gsub(" ", "_", db1$scientificName[i])
  
  # Imprimir el nombre de la especie para depuración
  print(paste("Checking species:", species))
  
  if(point %ni% c("BEP1 ", paste0("BF", 1:6), paste0("BP", 1:6), NA, "sep-01", "sep-02", "sep-03",
                  "TAP1",
                  "FOREST_NEAR1", "FOREST_NEAR1A","FOREST_NEAR2", "FOREST_NEAR2A", "FOREST_NEAR3") &
     species %ni% c("", "Uroxys_sp._3", "Coprophanaeus_gr_pluto", "Canthon_colombianus",
                    "Ateuchus_sp._01H", "Canthon_fulgidus", "Ateuchus_sp._09H", "Scybalocanthon_arnaudi_",
                    "Dichotomius_boreus", "Ateuchus_sp._2", "Ateuchus_sp._08H", "Ateuchus_sp._10H",
                    "Gromphas_lemoinei_", "Onthophagus_sp._2J", "Onthophagus_landolti_", "Oxysternon_silenus_",
                    "non_spp_record", "Anisocanthon_villosus_", "Sylvicanthon_sp._03H", "Ateuchus_sp._05H",
                    "Onthophagus_curivicornis", "Onthophagus_sp.08H", "Canthidium_sp._4FE")) {
    
    if(point %ni% points_sf$point) { stop("point issue") }
    
    # Comprobar la existencia de la especie
    if(species %ni% dist$scientificName) { 
      stop(paste("name issue: species", species, "not found in dist"))
    }
    
    distRegion <- dist$Region[dist$scientificName == species]
    pointsRegion <- points_sf$dist_region[points$point == point]
    
    if(!grepl(pointsRegion, distRegion)) {
      counter1 <- counter1 + 1
      region_mismatch[counter1] <- i
    }
    
    distSlope <- dist$Slope[dist$scientificName == species]
    pointsSlope <- points_sf$dist_slope[points$point == point]
    
    if(grepl(pointsRegion, "Andean")) {
      pointsSlope <- gsub(" ", "", pointsSlope)
      distSlope <- gsub(" ", "", distSlope)
      pointsSlope <- unlist(strsplit(pointsSlope, ":"))
      distSlope <- unlist(strsplit(distSlope, "(:|;)"))
      
      if(pointsSlope[1] %ni% distSlope) {
        counter2 <- counter2 + 1
        slope_mismatch[counter2] <- i
      } else {
        n <- which(distSlope == pointsSlope[1])
        l <- length(distSlope)
        
        if(l < (n + 2)) {
          if(!grepl(pointsSlope[2], distSlope[l])) {
            counter2 <- counter2 + 1
            slope_mismatch[counter2] <- i
          }
        } else if (n + 2 <= l) {
          # Verifica en el rango de n+1 a n+2
          if(!any(grepl(pointsSlope[2], distSlope[(n + 1):(n + 2)]))) {
            counter2 <- counter2 + 1
            slope_mismatch[counter2] <- i
          }
        } else {
          # Manejo del caso donde l es menor que n+2 pero mayor que n
          if(!grepl(pointsSlope[2], distSlope[(n + 1)])) {
            counter2 <- counter2 + 1
            slope_mismatch[counter2] <- i
          }
        }
      }
    }
  }
}
region_mismatch_df <- data.frame(species = db1$scientificName[region_mismatch], 
                                 point = db1$point[region_mismatch],
                                 species_dist = NA, point_dist = NA)
for(i in 1:nrow(region_mismatch_df)){
  region_mismatch_df$species_dist[i] <- dist$Region[dist$scientificName == gsub(" ", "_", region_mismatch_df$species[i])]
  region_mismatch_df$point_dist[i] <- points_sf$dist_region[points_sf$point == region_mismatch_df$point[i]]
}


slope_mismatch_df <- data.frame(species = db1$scientificName[slope_mismatch], 
                                point = db1$point[slope_mismatch],
                                species_dist = NA, point_dist = NA)
for(i in 1:nrow(slope_mismatch_df)){
  slope_mismatch_df$species_dist[i] <- dist$Slope[dist$scientificName == gsub(" ", "_", slope_mismatch_df$species[i])]
  slope_mismatch_df$point_dist[i] <- points_sf$dist_slope[points_sf$point == slope_mismatch_df$point[i]]
}
slope_mismatch_df <- slope_mismatch_df[order(slope_mismatch_df$point_dist),]

mismatch_species <- unique(slope_mismatch_df$species)
