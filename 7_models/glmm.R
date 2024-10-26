# 
# Generalized Linear Mixed Model (GLMM) with glmer
    db5 <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db5_traits.RDS")
    db5$elev_standard_squared <- db5$elev_standard^2
    db5$subregion_species <- paste0(db5$subregion, "__", db5$scientificName)
    db5$cluster_species <- paste0(db5$cluster, "__", db5$scientificName)
# 
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
      abundance ~ pasture + combinations + elev_standard + elev_standard_squared +  
                  nest_guild + diet_range + activity + bodysize + legratio +  
                  (1 + pasture + elev_standard + elev_standard_squared + nest_guild + 
                  diet_range + activity + bodysize + legratio| scientificName) +
                  (1 | subregion_species) + 
                  (1 | cluster_species),
                  data = db5,
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
