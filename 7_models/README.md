# Hierarchical modeling of dung beetle abundance


This script [glmm_brms.R](./glmm_brms.R) implements the statistical models used to evaluate drivers of dung beetle abundance across sampling sites in Colombia.

The analyses are based on the unified dataset produced in 6_DB_dataset, which integrates dung beetle abundance records with environmental predictors, geographic constraints, and functional trait information derived from previous modules of the repository.

## Generalized linear mixed models (GLMM)

A series of generalized linear mixed models were first fitted using the glmmTMB package to explore alternative random-effects structures and assess model convergence under different hierarchical specifications.

These exploratory GLMMs evaluated models with:

- species-specific random intercepts
- species-specific random slopes for environmental predictors
- nested spatial random effects associated with subregions and sampling clusters
- interactions between land use (pasture) and functional traits

These models were used exclusively to guide the specification of the hierarchical structure and to identify model configurations that produced stable estimates.

# Bayesian hierarchical model

The final statistical analysis is implemented using a Bayesian hierarchical negative binomial model fitted with the brms package (Stan backend). This model estimates species-specific responses to environmental gradients and land-use change while accounting for hierarchical structure across species and spatial sampling units.

The Bayesian model corresponds to the statistical framework described in “*Bayesian species-specific abundance modeling and prediction*” and “*Full Bayesian hierarchical negative binomial mixed-effects model specification and validation*” in Supplementary Online Material (SOM).

## Model predictors

The environmental predictors and trait covariates used in the model are derived from previous modules of the repository:

- geographic range constraints and spatial regionalization [3_Geographic_range](../3_Geographic_range)
- functional trait data (morphometric and behavioral traits) [4_Traits](../4_Traits)
- standardized elevation gradient derived from species elevation ranges [5_Standardized_elevation](../5_Standardized_elevation)
- unified abundance dataset and derived covariates [6_DB_dataset](../6_DB_dataset)
