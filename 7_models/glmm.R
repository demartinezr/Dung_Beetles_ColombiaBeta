# Generalized linear mixed models and bayesian analysis
library(glmmTMB)
library(dplyr)
library(ggplot2)
# dataset with abundance by speices/point/day -> 923315 obs, 958 points, 243 species
    db5 <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db5_distance.RDS")
    db5$elev_standard_squared <- db5$elev_standard^2
    db5$subregion_species <- paste0(db5$subregion, "__", db5$scientificName)
    db5$cluster_species <- paste0(db5$cluster, "__", db5$scientificName)
    db5$distance_from_range_scaled2 <- as.vector(db5$distance_from_range_scaled2)
  # sum abundance day by species/point. Remove day covariate -> 241948 obs, 958 points, 243 species
    db6 <- db5 %>%  
      group_by(across(-c(day, p_d))) %>%
      summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")  
# remove palm and young forest 219452 obs, 870 points, 243 species
    db6 <- db6[db6$habitat !="Sy",]
    db6 <- db6[db6$habitat !="PALM",]
#
# Generalized Linear Mixed Model (GLMM) with glmmTMB. Abundance by species/point
#
# protocol based on https://fhernanb.github.io/libro_modelos_mixtos/pac-glmmTMB.html
# https://stackoverflow.com/questions/73135114/convergence-criteria-in-glmmtmb-what-are-my-options
#
# Model with species-specific random slopes and intercepts for all covariates
    db_mod2_glmmTMB <- glmmTMB(
      abundance ~ pasture + distance_from_range_scaled2 + elev_standard + elev_standard_squared +  
                  nest_guild + diet_range + activity + bodysize + legratio +  
                  (1 + pasture + elev_standard + elev_standard_squared + nest_guild + 
                     diet_range + activity + bodysize + legratio| scientificName) +
                  (1 | subregion_species) + 
                  (1 | cluster_species),
                  data = db6,
                  family = nbinom2,
                  control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))
                  )
    # summary of adjust model no convergence
    saveRDS(db_mod2_glmmTMB, "C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod2_glmmTMB.rds")
    summary(db_mod2_glmmTMB)
#   result mod2
#    Warning message:
#      In finalizeTMB(TMBStruc, obj, fit, h, data.tmb.old) :
#      Model convergence problem; non-positive-definite Hessian matrix. See vignette('troubleshooting')
#    Warning message:
#      In sqrt(diag(vcovs)) : Se han producido NaNs
#
#   Model with independent random intercepts for species, subregion, and species cluster
    db_mod3_glmmTMB <- glmmTMB(
      abundance ~ pasture + distance_from_range_scaled2 + elev_standard + elev_standard_squared +  
        nest_guild + diet_range + activity + bodysize + legratio +  
        (1 | scientificName) + 
        (1 | subregion_species) + 
        (1 | cluster_species),
      data = db6,
      family = nbinom2,
      control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))
    )
    # summary of adjust model ## Get convergence
    saveRDS(db_mod3_glmmTMB, "C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod3_glmmTMB.rds")
    summary(db_mod3_glmmTMB)
    #
#    
# Model with species-specific random intercepts and random slopes for pasture, 
# including random intercepts for subregion and species cluster     
    db_mod4_glmmTMB <- glmmTMB(
      abundance ~ pasture + distance_from_range_scaled2 + elev_standard + elev_standard_squared +  
        nest_guild + diet_range + activity + bodysize + legratio +  
        (1 + pasture | scientificName) + # 1 + elev_standard + elev_standard_squared + nest_guild + diet_range + activity + bodysize + legratio| scientificName
        (1 | subregion_species) + 
        (1 | cluster_species),
      data = db6,
      family = nbinom2,
      control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))
    )
    # summary of adjust model ## Get convergence
    saveRDS(db_mod4_glmmTMB, "C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod4_glmmTMB.rds")
    summary(db_mod4_glmmTMB)
  #
  #
  #
    db_mod5_glmmTMB <- glmmTMB(
      abundance ~ pasture + distance_from_range_scaled2 + elev_standard + elev_standard_squared +  
        nest_guild + diet_range + activity + bodysize + legratio +  
        (1 + pasture + elev_standard_squared | scientificName) + # 1 + elev_standard + nest_guild + diet_range + activity + bodysize + legratio| scientificName
        (1 | subregion_species) + 
        (1 | cluster_species),
      data = db6,
      family = nbinom2,
      control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))
    )
    # summary of adjust model ## Get convergence
    saveRDS(db_mod5_glmmTMB, "C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod5_glmmTMB.rds")
    summary(db_mod5_glmmTMB)
    #
    #
#
# Model with species-specific random intercepts and random slopes for pasture and elevation terms, 
# including random intercepts for subregion and species cluster
    db_mod6_glmmTMB <- glmmTMB(
      abundance ~ pasture + distance_from_range_scaled2 + elev_standard + elev_standard_squared +  
        nest_guild + diet_range + activity + bodysize + legratio +  
        (1 + pasture + elev_standard + elev_standard_squared | scientificName) + # 1 + elev_standard + nest_guild + diet_range + activity + bodysize + legratio| scientificName
        (1 | subregion_species) + 
        (1 | cluster_species),
      data = db6,
      family = nbinom2,
      control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))
    )
    # summary of adjust model ## Get convergence
    saveRDS(db_mod6_glmmTMB, "C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod6_glmmTMB.rds")
    summary(db_mod6_glmmTMB)
    #
# Model with species-specific random intercepts and slopes, including traits (nest_guild)
# in random effects (non-convergent)
    db_mod7_glmmTMB <- glmmTMB(
      abundance ~ pasture + distance_from_range_scaled2 + elev_standard + elev_standard_squared +  
        nest_guild + diet_range + activity + bodysize + legratio +  
        (1 + pasture + elev_standard + elev_standard_squared + nest_guild | scientificName) + # 1 + elev_standard + nest_guild + diet_range + activity + bodysize + legratio| scientificName
        (1 | subregion_species) + 
        (1 | cluster_species),
      data = db6,
      family = nbinom2,
      control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))
    )
    # summary of adjust model ## No convergence
    saveRDS(db_mod7_glmmTMB, "C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod7_glmmTMB.rds")
    summary(db_mod7_glmmTMB)
#    
# Model with species-specific random intercepts and slopes, including traits (bodysize)
# in random effects (non-convergent)
      db_mod8_glmmTMB <- glmmTMB(
      abundance ~ pasture + distance_from_range_scaled2 + elev_standard + elev_standard_squared +  
        nest_guild + diet_range + activity + bodysize + legratio +  
        (1 + pasture + elev_standard + elev_standard_squared + bodysize | scientificName) + # 1 + elev_standard + nest_guild + diet_range + activity + bodysize + legratio| scientificName
        (1 | subregion_species) + 
        (1 | cluster_species),
      data = db6,
      family = nbinom2,
      control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))
    )
    # summary of adjust model ## No convergence
    saveRDS(db_mod8_glmmTMB, "C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod7_glmmTMB.rds")
    summary(db_mod8_glmmTMB)
    
# Remove the improbable combinations for point and species range, 865 points, 243 species
    db7 <- db6 |>
      filter(combinations ==1)
#
# Model included a trait-pasture interactions
    db_mod9_glmmTMB <- glmmTMB(
      abundance ~ pasture + elev_standard + elev_standard_squared + # biophysical
        nest_guild + diet_range + activity + bodysize + legratio + # traits
        pasture * nest_guild + pasture * diet_range + pasture * activity + # trait-pasture interactions
        pasture * bodysize + pasture * legratio +
        (1 + pasture + elev_standard + elev_standard_squared | scientificName) + 
        (1 | subregion_species) + 
        (1 | cluster_species),
      data = db7,
      family = nbinom2,
      control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))
    )
saveRDS(db_mod9_glmmTMB, "C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod9_glmmTMB.rds")
summary(db_mod9_glmmTMB)

library(brms)
db_mod1 <-
  brm(abundance ~
        pasture + elev_standard + elev_standard_squared +  # biophysical
        nest_guild + diet_range + activity + bodysize + legratio + # traits
        pasture * nest_guild + pasture * diet_range + pasture * activity + # trait-pasture interactions
        pasture * bodysize + pasture * legratio +
        (1 + pasture + elev_standard + elev_standard_squared | scientificName) +
        (1 | subregion_species) + (1 | cluster_species),
      family = "negbinomial", data = db7,
      chains = 3, cores = 3, backend = 'cmdstanr',
      refresh = 10)

db_mod2 <- brm(
  abundance ~ pasture + elev_standard + elev_standard_squared +
    nest_guild + diet_range + activity + bodysize + legratio +
    pasture * nest_guild + pasture * diet_range + pasture * activity +
    pasture * bodysize + pasture * legratio +
    (1 + pasture + elev_standard + elev_standard_squared | scientificName) +
    (1 | subregion_species) + (1 | cluster_species),
  family = "negbinomial", data = db7,
  chains = 3, cores = 3, backend = 'rstan',
  refresh = 10
)
