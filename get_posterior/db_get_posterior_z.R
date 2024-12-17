library(brms)
library(dplyr)
library(tidyr)

db_mod_abundance <-readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod_abundance.rds")
db_mod <- as.data.frame(db_mod_abundance)

# Functions to extract posterior iteration for several quantities of interest
get_psi_components <- function(fit) { 
  linear_pred <- posterior_predict(fit, re_formula = NULL) 
  psi_components <- linear_pred[, 1] 
  return(psi_components)
}
get_theta_component <- function(fit) { 
  linear_pred <- posterior_predict(fit, re_formula = NULL) 
  theta_component <- linear_pred[, 2] 
  return(theta_component)
}
get_z_components <- function(fit) { 
  psi_comp <- get_psi_components(fit) 
  theta_comp <- get_theta_component(fit) 
  return(list(psi = psi_comp, theta = theta_comp))
}
get_prediction_components <- function(fit) { 
  psi_comp <- get_psi_components(fit) 
  theta_comp <- get_theta_component(fit) 
  coefficients <- fixef(fit) 
  hyperparameters <- ranef(fit) 
  return(list(psi = psi_comp, theta = theta_comp, coefficients = coefficients, hyperparameters = hyperparameters))
}
get_Z_probs <- function(fit) { 
  psi_comp <- get_psi_components(fit) 
  z_probs <- posterior_predict(fit) 
  return(list(psi = psi_comp, Z = z_probs))
}
get_data_rep <- function(fit) { 
  posterior_data <- posterior_predict(fit) 
  return(posterior_data)
}

# Extract and save quantities of interest from the fitted model `db_mod_abundance`
psi_components <- get_psi_components(db_mod_abundance) 
theta_component <- get_theta_component(db_mod_abundance) 
z_components <- get_z_components(db_mod_abundance) 
prediction_components <- get_prediction_components(db_mod_abundance) 
z_probs <- get_Z_probs(db_mod_abundance) 

data_rep <- get_data_rep(db_mod_abundance)
posterior_samples <- as.data.frame(fitted(db_mod_abundance))  
num_samples <- nrow(posterior_samples) 
num_rows_db7 <- nrow(db7)

# Combine the extracted components into a dataframe
# Assuming that species names, sampling points, and observed Q values are available in the model's data (db7)
output <- data.frame(
  fitted_psi = psi_components,  # Fitted psi components
  fitted_conditional_prob = prediction_components$theta,  # Conditional probability of detection
  posterior_prob_Z = z_probs$Z,  # Posterior probability for Z = 1
  conditional_detection_probs = theta_component  # Detection probability per visit
)

saveRDS(output, "C:/Users/Dell-PC/Dropbox/CO_DBdata/Analysis/posterior_iteration.rds")
################################################################################
# total adapted
#
# Funciones para extraer iteraciones posteriores para varias cantidades de interés para un modelo de abundancia

# Extraer los componentes de abundancia (predicciones posteriores)
# Obtener las predicciones posteriores del modelo, sin los efectos aleatorios
get_abundance_components <- function(fit) {
  posterior_pred <- posterior_predict(fit, re_formula = NULL)
  abundance_components <- posterior_pred
  return(abundance_components)
}

# Calcular la abundancia media posterior
get_abundance_mean <- function(fit) {
  abundance_components <- get_abundance_components(fit)
  mean_abundance <- apply(abundance_components, 2, mean)  # Promedio de la abundancia posterior por parámetro
  return(mean_abundance)
}

# Calcular intervalos de credibilidad para la abundancia (por ejemplo, al 95%)
get_abundance_interval <- function(fit, prob = 0.95) {
  abundance_components <- get_abundance_components(fit)
  lower_bound <- apply(abundance_components, 2, function(x) quantile(x, probs = (1 - prob) / 2))
  upper_bound <- apply(abundance_components, 2, function(x) quantile(x, probs = 1 - (1 - prob) / 2))
  abundance_interval <- rbind(lower_bound, upper_bound)  # Combina límites en una matriz
  return(abundance_interval)
}

# Calcular probabilidades de interés para la abundancia (por ejemplo, probabilidad de ser mayor que 0)
get_abundance_probabilities <- function(fit) {
  abundance_components <- get_abundance_components(fit)
  prob_greater_than_zero <- apply(abundance_components, 2, function(x) mean(x > 0))
  return(prob_greater_than_zero)
}

# Extraer y guardar todas las cantidades de interés
get_abundance_metrics <- function(fit) {
  abundance_components <- get_abundance_components(fit)
  mean_abundance <- get_abundance_mean(fit)
  abundance_interval <- get_abundance_interval(fit)
  prob_greater_than_zero <- get_abundance_probabilities(fit)
  
  # Unir los resultados en una lista para mayor facilidad de acceso
  return(list(
    abundance = abundance_components,  # Componentes de abundancia posterior
    mean_abundance = mean_abundance,  # Abundancia media posterior
    abundance_interval = abundance_interval,  # Intervalo de credibilidad de la abundancia
    prob_greater_than_zero = prob_greater_than_zero  # Probabilidad de que la abundancia sea mayor que 0
  ))
}

# Extraer y guardar las métricas de interés del modelo ajustado `db_mod_abundance`
abundance_metrics <- get_abundance_metrics(db_mod_abundance)

# Combinar las métricas en un dataframe para facilitar la interpretación
output <- data.frame(
  abundance = apply(abundance_metrics$abundance, 2, mean),  # Abundancia posterior promedio
  mean_abundance = abundance_metrics$mean_abundance,  # Abundancia media
  abundance_interval_lower = abundance_metrics$abundance_interval[1,],  # Límite inferior del intervalo
  abundance_interval_upper = abundance_metrics$abundance_interval[2,],  # Límite superior del intervalo
  prob_greater_than_zero = abundance_metrics$prob_greater_than_zero  # Probabilidad de que la abundancia sea mayor que 0
)
# Ver resultados
print(output)

###############################################################################
#
# Paquetes necesarios
library(brms)
library(tidyverse)

# Función para extraer componentes del modelo
get_abundance_components <- function(model) {
  # Extrae los efectos fijos y aleatorios del modelo
  fixed_effects <- as_draws_df(model, variable = "b_") %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")
  
  random_effects <- as_draws_df(model, variable = "r_") %>%
    pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")
  
  list(fixed_effects = fixed_effects, random_effects = random_effects)
}

# Función para calcular probabilidades posteriores para los puntos
get_posterior_abundance <- function(model, newdata) {
  posterior_predictions <- posterior_epred(model, newdata = newdata)
  posterior_summary <- posterior_summary(posterior_predictions)
  list(posterior_predictions = posterior_predictions, posterior_summary = posterior_summary)
}

# Función para obtener componentes de predicción
get_prediction_components <- function(model, newdata) {
  abundance_components <- get_abundance_components(model)
  
  predictions <- posterior_epred(model, newdata = newdata)
  prediction_summary <- posterior_summary(predictions)
  
  list(
    fixed_effects = abundance_components$fixed_effects,
    random_effects = abundance_components$random_effects,
    predictions = predictions,
    prediction_summary = prediction_summary
  )
}

# Función para simular un conjunto de datos posterior
get_data_rep <- function(model, newdata) {
  simulated_data <- posterior_predict(model, newdata = newdata)
  data_summary <- posterior_summary(simulated_data)
  list(simulated_data = simulated_data, data_summary = data_summary)
}

# Ejemplo de uso con el modelo db_mod_abundance
# Cambia `newdata` para representar los puntos de interés
newdata <- data.frame(
  pasture = c(0, 1), 
  elev_standard = seq(-2, 2, length.out = 10),
  elev_standard_squared = seq(-2, 2, length.out = 10)^2,
  nest_guild = "Paracoprid",
  diet_range = "coprophagous",
  activity = "diurnal",
  bodysize = 5,
  legratio = 0.3,
  scientificName = "Species1",
  subregion_species = "Subregion1",
  cluster_species = "Cluster1"
)


