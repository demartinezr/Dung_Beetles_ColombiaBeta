setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata/SIG/inputs")

##### Ingest and clean beetle data #####
#
library(readxl)
# Dung beetles data set project all projects Western cordillera, llanos, Paramos, beta diversity
  db <- read_excel("C:/Users/Dell-PC/Dropbox/CO_dbdata/abundance/Scarabaeinae_database_2024.xlsx", sheet = "Scarabaeinae_database_2024")
  db$scientificName1 <- gsub("_", " ", db$scientificName)
  db$decimalLongitude <- as.numeric(db$lon_all_points)
  db$decimalLatitude <- as.numeric(db$lat_all_points)
  db$coordinates <- paste(db$lon_all_points, db$lat_all_points, sep = "_")
  db2 <- db[, c("point", "day", "scientificName", "abundance")]
#
# All sites visited in the projects: Western cordillera, llanos, Paramos, beta diversity
  j_points <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/all_pts.rds")
  pts_add <- data.frame(point = "TAP1", lat = 5.997912968, lon = -73.220879,
                        site = "El Taladro", cluster = "cluster_TAP_1", birds = 0, beetles = 1,
                        habitat = "Pasture", natural = 0, paramo = 0, pasture = 1,other = 0,mixed_cluster = 0,
                        elev_ALOS = 2350, posix1 = NA, posix2 = NA, posix3 = NA, posix4 = NA, obs1 = NA,
                        obs2 = NA, obs3 = NA, obs4 = NA, oday1 = NA, oday2 = NA, oday3 = NA, oday4 = NA,
                        hps1 = NA, hps2 = NA, hps3 = NA, hps4 = NA, nv = 4, subregion = "subregion_bogotacito",
                        stringsAsFactors = FALSE
  )
  j_points <- rbind(j_points, pts_add)
  j_points$habitat[is.na(j_points$habitat)] <- db$habitat_all_points[match(j_points$point[is.na(j_points$habitat)], db$point)]

### Identify and fix points in db2 that don't exist in j_points
  unique(db2$point[!(db2$point %in% j_points$point)])
# Naming errors
  db2$point[db2$point == "SEP_1"] <- "SEP1"
  db2$point[db2$point == "SEP_2"] <- "SEP2"
  db2$point[db2$point == "SEP_3"] <- "SEP3"
# Superfluous points: TAP1 (located too close to forest) and CHA7-9 (not properly collected)
  db2 <- db2[!(db2$point %in% c("TAP1", "CHA7", "CHA8", "CHA9")), ]
### Inspect and fix points in j_point that don't exist in db2
  missing_pts <- unique(j_points$point[!(j_points$point %in% db2$point) & (j_points$habitat != "PALM") & (j_points$habitat != "Sy")])
  missing_pts
# Confirm that nothing unexpected is happening here
  missing_pts2 <- j_points$point[!(j_points$point %in% db2$point) & (j_points$natural | j_points$pasture)]
  missing_pts[!(missing_pts %in% missing_pts2)]
### Remove points with insufficient collection day info
  no_day_points <- unique(db2$point[db2$day == ""])
  no_day_points
  db2 <- db2[!(db2$point %in% no_day_points), ]
### Get all point-days and species-point-days that exist in the data
  point_day <- paste(db2$point, db2$day, sep = "__")
  species_point_day <- paste(db2$scientificName, point_day, sep = "___")
# remove rows without a species (Important that we already generated point_day)
  db2 <- db2[db2$scientificName != "NA", ]
##### Zero-fill #####
  # Get all species; all point-days
  species <- unique(db2$scientificName)
  point_days <- unique(point_day)
  # Get all species-point-days that exist in db2
  species_point_day_db2 <- paste(db2$scientificName, db2$point, db2$day, sep = "__")
  # Assemble data-frame with all possible species-point-days
  all_spd <- data.frame(point = NA, day = NA, 
                        scientificName = rep(species, length(point_days)), 
                        abundance = 0, 
                        point_day = rep(point_days, each = length(species)))
  species_point_day_all <- paste(all_spd$scientificName, all_spd$point_day, sep = "__")
  first_cols <- do.call(rbind, strsplit(all_spd$point_day, "__"))
  all_spd$point <- first_cols[,1]
  all_spd$day <- first_cols[,2]
  
  all_spd <- all_spd[, c("point", "day", "scientificName", "abundance")]
  
  head(all_spd)
  head(db2)
  
  db2_additions <- all_spd[!(species_point_day_all %in% species_point_day_db2), ]
# Confirm that the dimension makes sense
  nrow(all_spd) == nrow(db2_additions) + nrow(db2)  
  db3 <- rbind(db2, db2_additions)
##### Add point data #####  
  j_points <- j_points[,c("point", "lat", "lon", "site", "cluster", 
                          "habitat", "natural", "paramo", "pasture", 
                          "other", "mixed_cluster", "elev_ALOS", "subregion")]
  
  source("D:/Doctorado/Tesis/Topographic_Col/unidades_topograficas.R")
  
  points_sf <- st_as_sf(j_points, coords = c("lon", "lat"), crs = 4326)
  j_points$pt_slope <- j_points$pt_region <- NA
  
  j_points$pt_region[st_intersects(points_sf, snsm, sparse = F)] <- "SNS Marta"
  j_points$pt_region[!st_intersects(points_sf, snsm, sparse = F) &
                       ((!st_intersects(points_sf, amazon_orinoco, sparse = F)) |
                          (points_sf$elev_ALOS > 500))] <- "Andean"
  j_points$pt_region[grepl("leguizamo", points_sf$subregion)] <- "Amazon"
  j_points$pt_region[grepl("chiribiquete", points_sf$subregion)] <- "Amazon"
  j_points$pt_region[grepl("guaviare", points_sf$subregion)] <- "Amazon"
  j_points$pt_region[grepl("llanos", points_sf$subregion)] <- "Llanos"
  
  j_points$pt_slope[st_intersects(points_sf, amazon_orinoco, sparse = F) |
                      st_intersects(points_sf, catatumbo, sparse = F)] <- "EC: Eastern"
  j_points$pt_slope[st_intersects(points_sf, magdalena_east, sparse = F)] <- "EC: Western"
  j_points$pt_slope[st_intersects(points_sf, magdalena_west, sparse = F)] <- "CC: Eastern"
  j_points$pt_slope[st_intersects(points_sf, cauca_east, sparse = F)] <- "CC: Western"
  j_points$pt_slope[st_intersects(points_sf, cauca_west, sparse = F)] <- "WC: Eastern"
  j_points$pt_slope[st_intersects(points_sf, pacific, sparse = F)] <- "WC: Western"
  j_points$pt_slope[j_points$pt_region == "SNS Marta"] <- "SNSM"
  
  db4 <- merge(db3, j_points, by.x = "point", by.y = "point", all = FALSE)
  
  ##### Format biogeography data #####
  biogeography <- read_excel("C:/Users/Dell-PC/Dropbox/CO_DBdata/DB_distri-20-07-21.xlsx", sheet="DB_distri-20-07-21")
  biogeography$Slope <- gsub(" ", "", biogeography$Slope)
  
  slope_split <- strsplit(biogeography$Slope, split = ";")
  unique(unlist(slope_split))
  
  biogeography2 <- data.frame(scientificName = biogeography$scientificName, 
                              sp_elev_lower = biogeography$Lower, sp_elev_upper = biogeography$Upper,
                              sp_region_amazon = NA, sp_region_llanos = NA,
                              sp_region_caribbean = NA, sp_region_snsm = NA,
                              sp_region_andes = NA, sp_region_eastern = NA,
                              sp_slope_ECe = NA, sp_slope_ECw = NA,
                              sp_slope_CCe = NA, sp_slope_CCw = NA,
                              sp_slope_WCe = NA, sp_slope_WCw = NA,
                              sp_slope_SNSM = NA)

    as.integer2 <- function (x) {
    if (length(x) == 0) {
      return(0)
    } else {
      return(max(as.integer(x)))
    }
  }
  
  for (i in 1:nrow(biogeography)) {
    biogeography2$sp_region_eastern[i] <- as.integer2(grepl("Eastern lowlands|Amazon|Llanos", biogeography$Region[i]))
    biogeography2$sp_region_amazon[i] <- as.integer2(grepl("Amazon", biogeography$Region[i]))
    biogeography2$sp_region_llanos[i] <- as.integer2(grepl("Llanos", biogeography$Region[i]))
    biogeography2$sp_region_caribbean[i] <- as.integer2(grepl("Caribbean", biogeography$Region[i]))
    biogeography2$sp_region_snsm[i] <- as.integer2(grepl("SNS Marta", biogeography$Region[i]))
    biogeography2$sp_region_andes[i] <- as.integer2(grepl("Andean", biogeography$Region[i]))
    
    slope_data <- slope_split[[i]]
    
    ec <- slope_data[grep("EC:", slope_data)]
    biogeography2$sp_slope_ECe[i] <- as.integer2(grepl("Eastern", ec))
    biogeography2$sp_slope_ECw[i] <- as.integer2(grepl("Western", ec))
    
    cc <- slope_data[grep("CC:", slope_data)]
    biogeography2$sp_slope_CCe[i] <- as.integer2(grepl("Eastern", cc))
    biogeography2$sp_slope_CCw[i] <- as.integer2(grepl("Western", cc))
    
    wc <- slope_data[grep("WC:", slope_data)]
    biogeography2$sp_slope_WCe[i] <- as.integer2(grepl("Eastern", wc))
    biogeography2$sp_slope_WCw[i] <- as.integer2(grepl("Western", wc))
    
    if(biogeography2$sp_region_eastern[i] == 1) {
      biogeography2$sp_slope_ECe[i] <- 1
    }
    
    # lowland species always cross valleys
    if ((biogeography2$sp_slope_ECw[i] != biogeography2$sp_slope_CCe[i]) & (biogeography2$sp_elev_lower[i] < 300)) {
      biogeography2$sp_slope_ECw[i] <- biogeography2$sp_slope_CCe[i] <- 1
    }
    
    if ((biogeography2$sp_slope_CCw[i] != biogeography2$sp_slope_WCe[i]) & (biogeography2$sp_elev_lower[i] < 300)) {
      biogeography2$sp_slope_CCw[i] <- biogeography2$sp_slope_WCe[i] <- 1
    }
    # 
    # # Caribbean species all potentially enter the Magdalena and reach the SNSM (max lower for these species is 102)
    # if (biogeography2$sp_region_caribbean[i] == 1) {
    #   biogeography2$sp_slope_ECw[i] <- biogeography2$sp_slope_CCe[i] <- biogeography2$sp_slope_SNSM[i] <- 1
    # }
    # 
    # # no gaps
    # ci <- which(colnames(biogeography2) == "sp_slope_ECe")
    # cf <- which(colnames(biogeography2) == "sp_slope_WCw")
    # the_numbers <- biogeography2[i, ci:cf]
    # if(sum(the_numbers) > 0){
    #   first_n <- min(which(the_numbers == 1))
    #   last_n <- max(which(the_numbers == 1))
    #   the_numbers[first_n:last_n] <- 1
    #   biogeography2[i, ci:cf] <- the_numbers
    # }
  }
  biogeography2$sp_slope_SNSM <- biogeography2$sp_region_snsm
  db5 <- merge(db4, biogeography2, by = "scientificName", all = TRUE)
# saveRDS(db5, "C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db_zerofill.RDS")
#
# add available species/point combinations in range column, based on shapefiles
#  library(doParallel)
#  library(foreach)
  library(sf)
  library(dplyr)
#
  db5 <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db_zerofill.RDS")
  shapefiles_combined <-  readRDS("C:/Users/Dell-PC/Documents/geographic_range.rds")
#  shapefiles_list <- list.files("D:/Doctorado/Tesis/GBIF/contorno", pattern = "*.shp", full.names = TRUE)
#  shapefiles_combined <- lapply(shapefiles_list, st_read) %>% do.call(rbind, .)
#  saveRDS(shapefiles_combined, "geographic_range.rds")
#
#
  df_occurrence <- db5[,c("scientificName","point","lon", "lat")]
  df_occurrence <- st_as_sf(df_occurrence, coords = c("lon", "lat"), crs = 4326)
  df_occurrence$combinations <- 0
  selected_shapefiles <- readRDS("C:/users/Dell-PC/documents/geographic_range.rds")
  for (i in 1:nrow(selected_shapefiles)) {
    # Filtrar por nombre científico del shapefile
    scientific_match <- selected_shapefiles$scientific[i]
    # Extraer el polígono correspondiente
    polygon <- selected_shapefiles[i, ]
    # Encontrar las filas de df_occurrence que coincidan con el nombre científico
    df_occurrence_match <- df_occurrence %>% 
      filter(scientificName == scientific_match)
    # Identificar los puntos dentro del polígono
    contained <- st_within(df_occurrence %>% st_as_sf(coords = c("lon", "lat"), crs = 4326), polygon, sparse = FALSE)
    df_occurrence$combinations[which(df_occurrence$scientificName == scientific_match & rowSums(contained) > 0)] <- 1
    rm(scientific_match, polygon, df_occurrence_match, contained)
    gc()
    }
  
  # Actualizar la columna 'combinations' en db5
  db5$combinations <- df_occurrence$combinations
  
  
# saveRDS(db5, "C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db_zerofillcomb.RDS")
  # expand 250 m, upper and lower boundaries 
  db5$sp_elev_lower2 <- pmin(db5$sp_elev_lower, (db5$sp_elev_lower + db5$sp_elev_upper)/2 - 250)
  db5$sp_elev_upper2 <- pmax(db5$sp_elev_upper, (db5$sp_elev_lower + db5$sp_elev_upper)/2 + 250)
  ##### Check biogeography #####
  db5$abundance <- as.numeric(db5$abundance)
  db_obs <- db5[db5$abundance > 0, ]
  # View(db_obs[db_obs$pt_region == "Andean" & db_obs$sp_region_andes == 0, ])
  # View(db_obs[db_obs$pt_region == "Amazon" & db_obs$sp_region_amazon == 0, ])
  # View(db_obs[db_obs$pt_region == "Llanos" & db_obs$sp_region_llanos == 0, ])
  # View(db_obs[db_obs$pt_region == "SNS Marta" & db_obs$sp_region_snsm == 0, ])
  # sum(db_obs$pt_slope == "SNSM" & db_obs$sp_slope_SNSM == 0)
  # View(db_obs[db_obs$pt_slope == "EC: Eastern" & db_obs$sp_slope_ECe == 0,])
  # View(db_obs[db_obs$pt_slope == "EC: Western" & db_obs$sp_slope_ECw == 0,])
  # sum(db_obs$pt_slope == "CC: Eastern" & db_obs$sp_slope_CCe == 0)
  # sum(db_obs$pt_slope == "CC: Western" & db_obs$sp_slope_CCw == 0)
  # sum(db_obs$pt_slope == "WC: Eastern" & db_obs$sp_slope_WCe == 0)
  # sum(db_obs$pt_slope == "WC: Western" & db_obs$sp_slope_WCw == 0)
  
##### Biogeographic clipping #####
  db5$elev_standard <- 2*(db5$elev_ALOS - db5$sp_elev_lower2)/(db5$sp_elev_upper2 - db5$sp_elev_lower2) - 1
  hist(db5$elev_standard[db5$abundance > 0])
  
  View(db5[db5$elev_standard > 3 & db5$abundance > 0, ])
#  
# add traits to dataset
  db_clipped <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db_clipped.RDS")
  beha <- read_excel("C:/Users/Dell-PC/Dropbox/CO_DBdata/traits/behaviour/DB_Distributions_traits_2024.xlsx", sheet="DB_Distributions_traits")
  beha <- beha[c("scientificName", "nest_guild", "diet_range", "activity")]
  morpho <- read.table("C:/Users/Dell-PC/Dropbox/CO_DBdata/traits/morphometrics/morphometrics_mean.txt")
  morpho <- morpho[morpho$scientificName %in% beha$scientificName, ]
  morphobeha <- merge(beha, morpho, by = "scientificName")
  db_clipped <- merge(db_clipped, morphobeha, by ="scientificName")
  db_clipped$nest_guild <- as.factor(db_clipped$nest_guild)
  db_clipped$diet_range <- as.factor(db_clipped$diet_range)
  db_clipped$activity <- as.factor(db_clipped$activity)
  db_clipped$subregion <- as.factor(db_clipped$subregion)
  db_clipped$cluster <- as.factor(db_clipped$cluster)
  db_clipped$bodysize <- as.numeric(db_clipped$bodysize)
  db_clipped$legratio <- as.numeric(db_clipped$legratio)
  saveRDS(db_clipped, "C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db_clipped_traits.RDS")  
  db_clipped <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db_clipped_traits.RDS")
#
#  add available/unavailable combinations (1/0) in "combined" column
  library(sf)
  library(dplyr)
  df_occurrence <- st_as_sf(db_clipped, coords = c("lon", "lat"), crs = 4326)
  # Create an empty data frame to store results
  results <- data.frame(
    point = df_occurrence$point,
    scientificName = df_occurrence$scientificName,
    combined = 0  # Inicializa todo en 0
  )
  shapefiles_list <- list.files("D:/Doctorado/Tesis/GBIF/contorno", pattern = "*.shp", full.names = TRUE)
  # Loop through each shapefile
  for (shapefile in shapefiles_list) {
    # Read the shapefile
    shapefile_species <- st_read(shapefile)
    # Check if occurrence points fall within the species range
    df_joined <- st_join(df_occurrence, shapefile_species, join = st_within)
     # Combine results
    results$combined[!is.na(df_joined$geometry)] <- 1
    # Liberar memoria
    rm(shapefile_species, df_joined)
    gc()
  }
  
  
#  
# Generalized Linear Mixed Model (GLMM) with glmer
  
    library(lme4)
  db_mod1_lme4 <- glmer(
    abundance ~ pasture + elev_standard + elev_standard_squared + nest_guild + 
      diet_range + activity + bodysize + legratio + subregion + cluster + 
      (1 + pasture + elev_standard + elev_standard_squared | scientificName) +
      (1 | subregion_species) + 
      (1 | cluster_species),
    data = db_clipped,
    family = negative.binomial(2),
    control = glmerControl(optCtrl = list(maxfun = 50000)),
    )
  # summary of adjust model
  summary(db_mod1_lme4)
  #
# Generalized Linear Mixed Model (GLMM) with glmmTMB
# protocol based on https://fhernanb.github.io/libro_modelos_mixtos/pac-glmmTMB.html
  library(glmmTMB)
  db_mod1_glmmTMB <- glmmTMB(
    abundance ~ pasture + elev_standard + elev_standard_squared + nest_guild + 
      diet_range + activity + bodysize + legratio +  
      (1 + pasture + elev_standard + elev_standard_squared | scientificName) +
      (1 | subregion_species) + 
      (1 | cluster_species),
    data = db_clipped,
    family = nbinom2,
    control=glmmTMBControl(optCtrl = list(maxit = 1000000)),
  )
  # summary of adjust model
  summary(db_mod1_glmmTMB)
  ##### brms model #####
  library(brms)
  db_clipped$elev_standard_squared <- db_clipped$elev_standard^2
  db_clipped$subregion_species <- paste0(db_clipped$subregion, "__", db_clipped$scientificName)
  db_clipped$cluster_species <- paste0(db_clipped$cluster, "__", db_clipped$scientificName)
  db_clipped <- na.omit(db_clipped)
  db_clipped$abundance <- as.numeric(db_clipped$abundance)
  
  db_mod1 <- 
    brm(abundance ~ pasture + elev_standard + elev_standard_squared + day + 
          (1 + pasture + elev_standard + elev_standard_squared + day | scientificName) +
          (1 | subregion_species) + (1 | cluster_species), 
        family = "negbinomial", data = db_clipped, 
        chains = 3, cores = 3, backend = 'cmdstanr',
        refresh = 10)
  