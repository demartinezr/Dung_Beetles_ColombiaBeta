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
get_prediction_components_abundance <- function(draws) {
  # Efectos fijos (coeficientes)
  fixed_effects <- draws[, grep("^b_", colnames(draws))]
  
  # Efectos aleatorios (desviaciones estándar y correlaciones)
  random_effects <- draws[, grep("^sd_", colnames(draws))]
  correlations <- draws[, grep("^cor_", colnames(draws))]
  
  # Interacciones (productos de coeficientes)
  interaction_terms <- draws[, grep(":", colnames(draws))]
  
  # Crear la lista de componentes
  components <- list(
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    correlations = correlations,
    interaction_terms = interaction_terms
  )
  
  return(components)
}

# Aplicar la función a los draws
prediction_components <- get_prediction_components_abundance(draws)

################################################################################
library(brms)
library(tidyr)
library(ggplot2)
library(ggokabeito)
library(ggthemes)

predicted <- posterior_predict(db_mod_abundance)
expected <- posterior_epred(db_mod_abundance)
linear_pred <- posterior_linpred(db_mod_abundance)
