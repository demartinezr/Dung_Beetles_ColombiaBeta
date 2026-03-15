# Predicting dung beetle abundance across Colombia

This script [db_prediction.R](./db_prediction.R) generates spatial predictions of dung beetle abundance across mainland Colombia using the Bayesian hierarchical model fitted in the repository.

Predictions are produced on a 2 × 2 km grid covering the country and estimate expected species abundance under alternative land-use conditions (forest vs pasture).

This prediction framework corresponds to the methods described in Methods "*Predicting species-specific abundance*", and "*Predictive modeling framework for dung beetle abundance*" in Supplementary Online Material

The script constructs the spatial prediction dataset and applies the Bayesian hierarchical model to estimate species-specific dung beetle abundance across Colombia.

Specifically, it:

- generates the national 2 × 2 km spatial grid
- integrates environmental predictors and species information
- restricts predictions to each species’ geographic distribution
- generates abundance predictions using posterior draws from the Bayesian model
- simulates abundance under forest and pasture land-use scenarios

These procedures produce spatially explicit predictions of dung beetle abundance while propagating parameter uncertainty from the Bayesian hierarchical model.
