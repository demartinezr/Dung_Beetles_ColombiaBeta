resultados <- read.table("./subregions_cluster.txt", header=TRUE)

library(brms)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

data_long <- pivot_longer(resultados, cols = c(clusters_natural, clusters_pasture), 
                          names_to = "habitat_type", values_to = "cluster_count")
data_long$subregion <- gsub("subregion_", "", data_long$subregion)

data_long$subregion <- factor(data_long$subregion, 
                              levels = data_long %>%
                                group_by(subregion) %>%
                                summarise(total_clusters = sum(cluster_count)) %>%
                                arrange(desc(total_clusters)) %>%
                                pull(subregion))

ggplot(data_long, aes(x = subregion, y = cluster_count, fill = habitat_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("clusters_natural" = "#4CAF50", "clusters_pasture" = "#FFC107"),
                    labels = c("Natural", "Pasture")) +
  labs(title = "Cluster distribution by habitat and subregion",
       x = "Subregion", y = "Number of clusters", fill = "Habitat type") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = c(4, 8, 12), color = "gray", alpha = 0.6, linetype = "dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
###############################################################################
db_mod1 <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod1.rds")
coef_modelo <- fixef(db_mod1) # coeficientes fijos
cred_intervals <- posterior_interval(db_mod1, prob = 0.95) # intervalos de credibilidad
rownames(cred_intervals) <- gsub("^b_", "", rownames(cred_intervals))
filas_comunes <- intersect(rownames(coef_modelo), rownames(cred_intervals))
coef_modelo <- coef_modelo[filas_comunes, , drop = FALSE]
cred_intervals <- cred_intervals[filas_comunes, , drop = FALSE]

coef_data <- data.frame(
  Variable = rownames(coef_modelo),
  Estimacion = coef_modelo[, "Estimate"],
  CI_lower = cred_intervals[, "2.5%"],
  CI_upper = cred_intervals[, "97.5%"]
)

coef_data$Significativo <- ifelse(coef_data$CI_lower > 0 | coef_data$CI_upper < 0, "Significant", "Not Significant")


ggplot(coef_data, aes(x = Estimacion, y = Variable, xmin = CI_lower, xmax = CI_upper, color = Significativo)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", alpha = 0.4, lwd=1) + # Línea en x = 0
    geom_point(size=3) + 
  geom_errorbarh(height = 0.6, linewidth = 0.8) +
  scale_color_manual(values = c("Significant" = "darkblue", "Not Significant" = "darkred")) +
  labs(title = "Fixed effects: coefficients and their credibility intervals",
    x = "Estimate", y = "Predictors") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 14), 
    axis.text = element_text(size = 14)
  )

################################################################################
library(brms)
library(tidyr)
library(ggplot2)
library(dplyr)
library(ggdist)

db_mod1 <- readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod1.rds")
# Extract random effects with their credibility intervals
random_effects <- ranef(db_mod1)

cluster_species <- as.data.frame(ranef(db_mod1)$cluster_species)
cluster_species$levels <- rownames(ranef(db_mod1)$cluster_species)
cluster_species_long <- cluster_species %>%
  pivot_longer(cols = -levels, 
               names_to = "parameter", 
               values_to = "values")
cluster_species_long$effect <- "cluster_species"

scientificName <- as.data.frame(ranef(db_mod1)$scientificName)
scientificName$levels <- rownames(ranef(db_mod1)$scientificName)
scientificName_long <- pivot_longer(
  scientificName,
  cols = -levels, #
  names_to = "parameter", 
  values_to = "values"
)
scientificName_long$effect <- "scientificName"

subregion_species <- as.data.frame(ranef(db_mod1)$subregion_species)
subregion_species$levels <- rownames(ranef(db_mod1)$subregion_species)
subregion_species_long <- pivot_longer(
  subregion_species,
  cols = -levels,   
  names_to = "parameter",  
  values_to = "values"  
)
subregion_species_long$effect <- "subregion_species"

random_effects_merge <- as.data.frame(bind_rows(cluster_species_long, 
                            scientificName_long, 
                            subregion_species_long))

# Filter data for intercept and specific parameters
filtered_data1 <- random_effects_merge %>%
  filter(parameter == "Estimate.Intercept")

filtered_data2 <- random_effects_merge %>%
  filter(parameter %in% c("Estimate.pasture1", "Estimate.elev_standard", "Estimate.elev_standard_squared"))

# Combine both filtered datasets into a single dataframe
combined_data <- bind_rows(filtered_data1, filtered_data2)

# Organize the Y-axis alphabetically
combined_data$effect_parameter <- factor(
  interaction(combined_data$effect, combined_data$parameter), 
  levels = sort(unique(interaction(combined_data$effect, combined_data$parameter)))
)

# Create ridgeline density plot for random effects, intercepts, and selected parameters
ggplot(combined_data, aes(x = values, y = effect_parameter, fill = effect_parameter)) +
  geom_density_ridges(alpha = 0.6, scale = 6, show.legend = FALSE) +  # Ridgeline densities
  scale_x_continuous(expand = c(0.1, 0)) +  # Avoid extra space on the X-axis
  theme_classic() +
  labs(
    title = "Ridgeline density Plot of estimated intercepts and specific parameters",
    x = "Estimates",
    y = "Random effects and parameters"
  ) +
  theme(
    panel.grid.major.y = element_blank(),  # Remove grid lines on the Y-axis
    axis.line.y = element_blank(),
    axis.text.y = element_text(size = 12),  # Improve Y-axis text readability
    plot.title = element_text(size = 16, hjust = 0.5),  # Center the title
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )
################################################################################
library(ggplot2)
library(reshape2)

m1 <- summary(db_mod1)
str(m1$fixed)

fixed_params <- data.frame(
  Parameter = rownames(m1$fixed),
  Rhat = m1$fixed$Rhat,
  Bulk_ESS = m1$fixed$Bulk_ESS,
  Tail_ESS = m1$fixed$Tail_ESS
) %>%
  pivot_longer(cols = c(Rhat, Bulk_ESS, Tail_ESS), 
               names_to = "Metric", 
               values_to = "Value")

# Create a faceted plot
ggplot(fixed_params, aes(x = Parameter, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Rhat" = "skyblue", "Bulk_ESS" = "orange", "Tail_ESS" = "darkgreen")) +
  labs(
    title = "Model diagnostics: Rhat, bulk ESS, and tail ESS",
    x = "Parameter",
    y = "Value"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 12),
    legend.position = "none"
  ) +
  # Add a reference line for Rhat (only for that facet)
  geom_hline(data = subset(fixed_params, Metric == "Rhat"), aes(yintercept = 1), 
             linetype = "dashed", color = "darkred", size = 0.8) +
  geom_hline(data = subset(fixed_params, Metric == "Bulk_ESS"), aes(yintercept = 400),
           linetype = "dotted", color = "darkred", size = 0.8) +
  geom_hline(data = subset(fixed_params, Metric == "Tail_ESS"), aes(yintercept = 400),
             linetype = "dotted", color = "darkred", size = 0.8)
################################################################################
db_mod1 <-readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod_abundance.rds")
# Marginal predictions using fixed effects
marginal_predictions <- fitted(db_mod1)

# Predictions with random effects
predictions_with_random_effects <- predict(db_mod1)

# Get the marginal effects
marginals <- conditional_effects(db_mod1)

# Plot the marginal effects for the "pasture" variable
plot(marginals, "pasture")

# Predictions with random effects (including random effects)
predictions_with_random_effects <- predict(db_mod1, re.form = NULL)

# Compare predictions
comparison_predictions <- data.frame(
  Observed = db7$abundance,  # Replace with your observed variable
  Fixed_Predictions = marginal_predictions[,1],  # Extract predictions from the correct column
  Random_Predictions = predictions_with_random_effects[,1]  # Extract predictions from the correct column
)

# Calculate prediction errors for both types of predictions
fixed_prediction_error <- comparison_predictions$Fixed_Predictions - comparison_predictions$Observed
random_prediction_error <- comparison_predictions$Random_Predictions - comparison_predictions$Observed

# Summary of the prediction errors
summary(fixed_prediction_error)
summary(random_prediction_error)

# Plot predictions with fixed and random effects
ggplot(comparison_predictions, aes(x = Observed)) +
  geom_point(aes(y = Fixed_Predictions), color = "darkblue", alpha = 0.5) +
  geom_point(aes(y = Random_Predictions), color = "darkred", alpha = 0.5) +
  labs(
    title = "Prediction comparison (Fixed effects vs Random effects)",
    x = "Observed abundance",
    y = "Prediction"
  ) +
  theme_classic()

# Add error columns to the comparison dataframe
comparison_predictions$error_fixed <- fixed_prediction_error
comparison_predictions$error_random <- random_prediction_error

# Plot prediction errors (fixed vs random)
ggplot(comparison_predictions, aes(x = Observed)) +
  geom_point(aes(y = error_fixed), color = "darkblue", alpha = 0.5) +
  geom_point(aes(y = error_random), color = "darkred", alpha = 0.5) +
  labs(
    title = "Prediction errors (Fixed vs Random)",
    x = "Observed abundance",
    y = "Prediction error"
  ) +
  theme_classic()
################################################################################
residuals <- residuals(db_mod1)
ggplot(data.frame(residuals), aes(x = Estimate)) +
  geom_histogram(binwidth = 0.1, fill = "lightblue", color = "black") +
  labs(title = "Distribución de los residuos", x = "Residuos", y = "Frecuencia")

# Calculate residuals (if you haven't done it yet)
residuals <- db7$abundance - marginal_predictions[,1]  # Replace 'db7$abundance' with your observed variable

# Create a dataframe with the predictions and residuals
data_plot <- data.frame(Predictions = marginal_predictions[,1], Residuals = residuals)

# Plot
ggplot(data_plot, aes(x = Predictions, y = Residuals)) +
  geom_point(color = "darkblue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Predictions", x = "Predictions", y = "Residuals") +
  theme_classic()

# Usar validación cruzada bayesiana con loo
library(loo)
loo_model <- loo(db_mod1)
loo_model

################################################################################
# 1. Prepare datasets for forest and pasture
forest_data <- db7 %>%
  mutate(pasture = 0, natural = 1)

pasture_data <- db7 %>%
  mutate(pasture = 1, natural = 0)

# 2. Generate predictions using posterior_predict()
library(brms)
forest_pred <- posterior_predict(db_mod1, newdata = forest_data)
pasture_pred <- posterior_predict(db_mod1, newdata = pasture_data)

# 3. Calculate sensitivity for each species
mean_abundance_forest <- apply(forest_pred, 2, mean)
mean_abundance_pasture <- apply(pasture_pred, 2, mean)
sensitivity <- log(mean_abundance_forest / mean_abundance_pasture)

species_names_db7 <- db7$scientificName
colnames(forest_pred) <- species_names_db7
# 4. Analyze sensitivity distribution
sensitivity_percentiles <- quantile(sensitivity, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
sensitivity_by_species <- colMeans(sensitivity, na.rm = TRUE) 

sensitivity_df <- data.frame(
  species = colnames(forest_pred),
  sensitivity = sensitivity)


ggplot(sensitivity_summary, aes(x = mean_sensitivity)) +
  geom_histogram(binwidth = 0.1, fill = "lightgray", color = "black") +
  geom_vline(xintercept = sensitivity_percentiles[1], col = "darkred", lwd = 1, lty = 2) +  # Percentil 25
  geom_vline(xintercept = sensitivity_percentiles[2], col = "darkgreen", lwd = 1, lty = 2) +  # Percentil 50
  geom_vline(xintercept = sensitivity_percentiles[3], col = "darkblue", lwd = 1, lty = 2) +   # Percentil 75
  labs(
    title = "Sensitivity distribution to forest conversion",
    x = "Sensitivity Log(forest abundance/pasture abundance)",
    y = "Frequency"
  ) +
  theme_classic() + 
  theme(legend.position = "none") +
  scale_y_continuous(expand = c(0, 0))


ggplot(sensitivity_summary, aes(x = reorder(species, mean_sensitivity), y = mean_sensitivity)) +
  geom_bar(stat = "identity", fill = "lightgray", color="black", size = 0.3) +  # Bordes en negro para mayor contraste
  coord_flip() +  # Rotar el eje para facilitar la lectura
  labs(
    title = "Mean sensitivity across species",
    x = "Species",
    y = "Mean Sensitivity (Log Scale)"
  ) +
  theme_classic() +  # Usa un tema limpio con texto base más grande
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),  # Título centrado y en negrita
    axis.text.y = element_blank(),  # Tamaño del texto de especies ajustado
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10),  # Tamaño del texto en el eje X ajustado
    axis.title.x = element_text(),  # Negrita en el título del eje X
    axis.title.y = element_text()  # Negrita en el título del eje Y
  ) +
  scale_y_continuous(expand = c(0, 0)) +  # Ajustar límites del eje Y
  geom_hline(yintercept = 0, color = "darkgray", linetype = "dashed") 

#######################################################################

min_sensitivity <- min(sensitivity[is.finite(sensitivity)], na.rm = TRUE)
max_sensitivity <- max(sensitivity[is.finite(sensitivity)], na.rm = TRUE)

# Calcular el promedio y manejar NA o Inf asignando mínimo o máximo
sensitivity_summary2 <- sensitivity_summary %>%
  mutate(
    mean_sensitivity = case_when(
      is.na(mean_sensitivity) ~ min_sensitivity,  # Reemplazar NA con el valor mínimo
      mean_sensitivity == Inf ~ max_sensitivity, # Reemplazar Inf con el valor máximo
      mean_sensitivity == -Inf ~ min_sensitivity, # Reemplazar -Inf con el valor mínimo
      TRUE ~ mean_sensitivity                     # Mantener valores válidos
    )
  )

# Estandarización de la sensibilidad a un rango de -1 a 1
sensitivity_summary2$scaled_sensitivity <- (sensitivity_summary2$mean_sensitivity - min(sensitivity_summary2$mean_sensitivity, na.rm = TRUE)) / 
  (max(sensitivity_summary2$mean_sensitivity, na.rm = TRUE) - min(sensitivity_summary2$mean_sensitivity, na.rm = TRUE)) * 2 - 1

# Graficar la distribución de la sensibilidad estandarizada
ggplot(sensitivity_summary2, aes(x = reorder(species, scaled_sensitivity), y = scaled_sensitivity)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +  # Rotar el gráfico para facilitar la lectura
  labs(
    title = "Scaled sensitivity across Species",
    x = "Species",
    y = "Scaled sensitivity"
  ) +
  theme_classic() +  # Usa un tema limpio con texto base más grande
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),  # Título centrado y en negrita
    axis.text.y = element_blank(),  # Tamaño del texto de especies ajustado
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10),  # Tamaño del texto en el eje X ajustado
    axis.title.x = element_text(),  # Negrita en el título del eje X
    axis.title.y = element_text()  # Negrita en el título del eje Y
  ) +
  scale_y_continuous(expand = c(0, 0)) +  # Ajustar límites del eje Y
  geom_hline(yintercept = 0, color = "darkgray", linetype = "dashed")   
