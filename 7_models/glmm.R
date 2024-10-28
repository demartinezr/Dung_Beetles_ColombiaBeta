# Generalized linear mixed models and bayesian analysis
library(glmmTMB)
library(dplyr)
library(ggplot2)
# dataset
    db5 <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db5_distance.RDS")
    db5$elev_standard_squared <- db5$elev_standard^2
    db5$subregion_species <- paste0(db5$subregion, "__", db5$scientificName)
    db5$cluster_species <- paste0(db5$cluster, "__", db5$scientificName)
    db5$distance_from_range_scaled2 <- as.vector(db5$distance_from_range_scaled2)
#
# Generalized Linear Mixed Model (GLMM) with glmmTMB
# protocol based on https://fhernanb.github.io/libro_modelos_mixtos/pac-glmmTMB.html
    db_mod1_glmmTMB <- glmmTMB(
      abundance ~ pasture + distance_from_range_scaled2 + elev_standard + elev_standard_squared +  
                  nest_guild + diet_range + activity + bodysize + legratio +  
                  (1 + pasture + elev_standard + elev_standard_squared + nest_guild + 
                  diet_range + activity + bodysize + legratio| scientificName) +
                  (1 | subregion_species) + 
                  (1 | cluster_species),
                  data = db5,
                  family = nbinom2,
                  control=glmmTMBControl(optCtrl = list(maxit = 1000000))
                  )
    # summary of adjust model
    saveRDS(db_mod1_glmmTMB, "C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod1_glmmTMB.rds")
    summary(db_mod1_glmmTMB)
#
##### brms model (J. Socolar)
library(brms)
db_mod1 <- 
  brm(abundance ~ pasture + combinations + elev_standard + elev_standard_squared +  
        nest_guild + diet_range + activity + bodysize + legratio +
        (1 + pasture + elev_standard + elev_standard_squared + nest_guild + 
           diet_range + activity + bodysize + legratio| scientificName) +
        (1 | subregion_species) + (1 | cluster_species), 
      family = "negbinomial", data = db5, 
      chains = 3, cores = 3, backend = 'cmdstanr',
      refresh = 10)
#
#
