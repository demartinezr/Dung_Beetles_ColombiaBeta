# ------------------------------------------------------------------------------
# Hierarchical modeling of dung beetle abundance
#
# This script fits generalized linear mixed models (GLMMs) and Bayesian
# hierarchical models to evaluate drivers of dung beetle abundance across
# sampling sites in Colombia.
#
# The dataset used here is produced in module `6_DB_dataset`, which integrates
# abundance data, environmental predictors, species geographic ranges,
# biogeographic constraints, and functional trait information.
#
# Several GLMMs are first fitted using the package `glmmTMB` to explore model
# structure and random-effects specifications. These exploratory models were
# used to evaluate convergence and determine an appropriate random-effects
# structure for the final hierarchical model.
#
# The final statistical analysis is implemented using a Bayesian hierarchical
# negative binomial model fitted with the package `brms` (Stan backend).
#
# The model estimates species-specific responses to:
#   • land-use change (pasture vs natural habitats)
#   • elevation gradients
#   • functional traits
#
# while accounting for hierarchical structure across:
#   • species
#   • sampling subregions
#   • spatial sampling clusters
# ------------------------------------------------------------------------------

setwd("C:/Users/PC/Dropbox/CO_DBdata")
library(glmmTMB)
library(dplyr)
library(ggplot2)
# ------------------------------------------------------------------------------
# Import the dataset produced in module `6_DB_dataset`.
# This dataset contains abundance observations linked to environmental
# predictors, geographic constraints, and functional trait information.
# speices/point/day -> 923315 obs, 958 points, 243 species
    db5 <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db5_distance.RDS")
    
# ------------------------------------------------------------------------------  
    # Create derived predictors and grouping variables used in the hierarchical
    # models (e.g., squared elevation and nested spatial grouping factors).    
    db5$elev_standard_squared <- db5$elev_standard^2
    db5$subregion_species <- paste0(db5$subregion, "__", db5$scientificName)
    db5$cluster_species <- paste0(db5$cluster, "__", db5$scientificName)
    db5$distance_from_range_scaled2 <- as.vector(db5$distance_from_range_scaled2)
  # sum abundance day by species/point. Remove day covariate -> 241948 obs, 958 points, 243 species
    db6 <- db5 %>%  
      group_by(across(-c(day, p_d))) %>%
      summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")  
  # Remove habitat types not included in the statistical analysis
  # (young secondary forest and oil palm plantations Gilroy et al. 2014, 2015).
  # remove palm and young forest 219452 obs, 870 points, 243 species
    db6 <- db6[db6$habitat !="Sy",]
    db6 <- db6[db6$habitat !="PALM",]
#
############################ Exploratory GLMM models ###########################
# Fit a series of generalized linear mixed models (GLMMs) using `glmmTMB`
# to explore alternative random-effects structures and evaluate model
# convergence prior to implementing the Bayesian hierarchical model.
#
# These models are not used for inference in the manuscript but helped
# determine the appropriate hierarchical structure for the final model.
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
    
######################## Geographic range filtering ############################
# Remove species–site combinations outside the estimated geographic
# ranges of species. Range constraints were derived in module
# `3_Geographic_range`. 865 points, 243 species
    db7 <- db6 |>
      filter(combinations ==1)
#
######################## Trait–environment interactions ########################
# Test models including interactions between land use (pasture) and
# functional traits to evaluate potential trait-mediated responses
# to land-use change.
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
#
######################## Bayesian hierarchical model ###########################
# Fit the final hierarchical model using `brms` (Stan backend).
#
# The model assumes a negative binomial distribution for abundance
# and estimates species-specific responses to environmental predictors
# and land-use change while accounting for spatial and taxonomic
# hierarchical structure.
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

#saveRDS(db_mod1, "C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod1.rds")

db_mod1 <-readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod1.rds")
summary(db_mod1)

######################## Final model used in the manuscript ####################
# This model corresponds to the Bayesian hierarchical model reported
# in the main manuscript and Supplementary Online Material.

class(db7$pasture)
db7$pasture <- as.factor(db7$pasture)
# saveRDS(db7, "C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/db7_available.RDS")
db_mod_abundance <-
  brm(abundance ~
        pasture + elev_standard + elev_standard_squared +  # biophysical
        nest_guild + diet_range + activity + bodysize + legratio + # traits
        pasture * nest_guild + pasture * diet_range + pasture * activity + # trait-pasture interactions
        pasture * bodysize + pasture * legratio +
        (1 + pasture + elev_standard + elev_standard_squared | scientificName) +
        (1 | subregion_species) + (1 | cluster_species),
      family = "negbinomial", data = db7,
      chains = 3, cores = 4, backend = 'cmdstanr',
      sample_file = "C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod_abundance.csv",
      refresh = 10)

#saveRDS(db_mod_abundance, "C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod_abundance.rds")
db_mod_abundance <-readRDS("C:/Users/PC/Dropbox/CO_DBdata/db_mod_abundance.rds")

summary(db_mod_abundance)
prior_summary(db_mod_abundance)