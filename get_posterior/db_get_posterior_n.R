library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)

db_mod_abundance <-readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod_abundance.rds")
draws <- as.data.frame(posterior::as_draws_df(db_mod_abundance))
n_info <- db_mod_abundance$data
iter <- draws$.iteration
# predicted <- posterior_predict(db_mod_abundance, re_formula = NULL)

# Function to get the species-specific components of the abundance model (log_lambda)
get_lambda_components <- function(draws, iter, n_info) {
  species <- n_info[!duplicated(n_info$scientificName),]
  # get the expectation at "pasture = 0" (base abundance when pasture is 0)
  log_lambda_0 <- 
    as.numeric(draws[iter, "b_Intercept"]) +
    as.numeric(draws[iter, paste0("r_scientificName[", species$scientificName, ",Intercept]")]) 
    # Get the effect of pasture on log_lambda (log-transformed abundance)
  log_lambda_pasture_offset <- 
    as.numeric(draws[iter, "b_pasture1"]) +
    as.numeric(draws[iter, paste0("r_scientificName[", species$scientificName, ",pasture1]")])
 
  return(data.frame(scientificName = species$scientificName, 
                    log_lambda_0 = log_lambda_0,
                    log_lambda_pasture_offset = log_lambda_pasture_offset))
}

lambda_components <- lapply(1:nrow(draws), function(i) get_lambda_components(draws, i, n_info))

saveRDS(lambda_components, "C:/Users/Dell-PC/Dropbox/CO_DBdata/Analysis/posterior_lambda.rds")
lambda <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/Analysis/get_posterior/posterior_lambda.rds")

################################################################################
library(brms)
library(tidyr)
library(ggplot2)
library(ggokabeito)
library(ggthemes)

predicted <- posterior_predict(db_mod_abundance)
expected <- posterior_epred(db_mod_abundance)
linear_pred <- posterior_linpred(db_mod_abundance)

################################################################################
# Function to get the species-specific components of the abundance model (log_lambda)
get_lambda_components <- function(draws, iter, n_info) {
  # Get unique species
  species <- n_info[!duplicated(n_info$scientificName),]
  
  # Get the expectation at "pasture = 0" (base abundance when pasture is 0)
  log_lambda_0 <- 
    as.numeric(draws[iter, "b_Intercept"]) +
    as.numeric(draws[iter, paste0("r_scientificName[", species$scientificName, ",Intercept]")]) 
  
  # Get the effect of pasture on log_lambda (log-transformed abundance)
  log_lambda_pasture_offset <- 
    as.numeric(draws[iter, "b_pasture1"]) +
    as.numeric(draws[iter, paste0("r_scientificName[", species$scientificName, ",pasture1]")])
  
  # Create an empty dataframe to store the results
  result_df <- data.frame(
    scientificName = species$scientificName,
    log_lambda_0 = log_lambda_0,
    log_lambda_pasture_offset = log_lambda_pasture_offset
  )
  
  # Add random effects for clusters (each cluster effect should be a column per species)
  for (species_name in species$scientificName) {
    # Find all clusters related to this species
    cluster_names <- unique(n_info$cluster_species[n_info$scientificName == species_name])
    
    for (cluster_name in cluster_names) {
      # Extract the random effect for this species and cluster
      cluster_effect_column <- paste0("cluster_", gsub("[^a-zA-Z0-9]", "_", cluster_name)) # Clean the name for valid column names
      cluster_effect_value <- as.numeric(draws[iter, paste0("r_cluster_species[", cluster_name, ",Intercept]")])
      
      # Add the cluster effect as a new column for the current species
      result_df[which(result_df$scientificName == species_name), cluster_effect_column] <- cluster_effect_value
    }
    
    # Find all subregions related to this species
    subregion_names <- unique(n_info$subregion_species[n_info$scientificName == species_name])
    
    for (subregion_name in subregion_names) {
      # Extract the random effect for this species and subregion
      subregion_effect_column <- paste0("subregion_", gsub("[^a-zA-Z0-9]", "_", subregion_name)) # Clean the name for valid column names
      subregion_effect_value <- as.numeric(draws[iter, paste0("r_subregion_species[", subregion_name, ",Intercept]")])
      
      # Add the subregion effect as a new column for the current species
      result_df[which(result_df$scientificName == species_name), subregion_effect_column] <- subregion_effect_value
    }
  }
  
  return(result_df)
}

lambda_components <- lapply(1:nrow(draws), function(i) get_lambda_components(draws, i, n_info))

saveRDS(lambda_components, "C:/Users/Dell-PC/Dropbox/CO_DBdata/Analysis/posterior_lambda_RES.rds")

