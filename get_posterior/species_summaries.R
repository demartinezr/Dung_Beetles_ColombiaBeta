# Fit the model using brms with the specified formula and options
# The model includes biophysical and trait variables, trait-pasture interactions, and random effects

library(brms)
library(posterior)

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
# Remove the improbable combinations for point and species range, 865 points, 243 species
db7 <- db6 |>
  filter(combinations ==1)

# Bayensian abundance model
db_mod_abundance <-readRDS("C:/Users/Dell-PC/Dropbox/CO_DBdata/db_mod_abundance.rds")

# Store posterior samples
posterior_samples <- as.array(db_mod_abundance)

# Convert posterior_samples into a data frame for easier manipulation
posterior_df <- as_draws_df(posterior_samples)

# Verify the variable names
variable_names <- dimnames(posterior_samples)$variable

# Extract species list and count
species_list <- unique(db7$scientificName)
species_count <- length(species_list)

# Initialize posterior estimates array
posterior_estimates <- matrix(NA, nrow = species_count, ncol = nrow(posterior_df))

# Iterate over species to extract random effects
for (sp_idx in seq_along(species_list)) {
  sp <- species_list[sp_idx]  # Current species
  
  # Generate column name for the random intercept
  random_effect_param <- paste0("r_scientificName[", sp, ",Intercept]")
  
  # Check if the column exists in posterior_df
  if (random_effect_param %in% colnames(posterior_df)) {
    # Extract random intercept for the species
    random_intercept <- posterior_df[[random_effect_param]]
    
    # Calculate posterior estimates for pasture effect
    posterior_estimates[sp_idx, ] <- random_intercept + posterior_df[["b_pasture1"]]
  } else {
    warning(paste("Random effect for species", sp, "not found in posterior samples."))
  }
}

# Name the rows of posterior estimates
rownames(posterior_estimates) <- species_list
saveRDS(posterior_estimates, "C:/Users/Dell-PC/Dropbox/CO_DBdata/Analysis/get_posterior/species_summaries.rds")

# Extract random effects coefficients from the model
fixef_summary <- summary(db_mod_abundance)$fixed
print(fixef_summary)

# Get posterior samples of the parameters
ranef_summary <- summary(db_mod_abundance)$random
print(ranef_summary)

# Obtener las muestras posterior de los parámetros
posterior_samples <- as_draws_df(db_mod_abundance)

# Analyze the effects of "pasture" and its interactions with other variables
pasture_effects <- posterior_samples %>%
  select(starts_with("b_pasture")) %>% # Select all parameters starting with "b_pasture"
  reframe(across(everything(), ~ quantile(., probs = c(0.05, 0.95, 0.1, 0.9, 0.025, 0.975, 0.5)))) %>% # Calculate quantiles
  t() %>% # Transpose the data frame for easier interpretation
  as.data.frame() # Convert to data frame

# Add meaningful column names to the quantile results
colnames(pasture_effects) <- c("q05", "q95", "q10", "q90", "q025", "q975", "mean")
rownames(pasture_effects) <- gsub("b_", "", colnames(select(posterior_samples, starts_with("b_pasture"))))

print("Efectos de 'pasture':")
print(pasture_effects)

# Identify species with significant effects of 'pasture'
# (those where the 95% credible interval does not cross zero)
pasture_significant_species <- pasture_effects %>%
  filter(q05 > 0 | q95 < 0)
# Print species with significant effects
print(pasture_significant_species)

# Print species with significant effects
hist(pasture_effects$mean, main = "Distribution of aAverage Pasture effects",
     xlab = "Average effect", col = "skyblue", border = "white")

################################################################################
posterior_predict_db <- posterior_predict(db_mod_abundance)
