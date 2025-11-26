This repository contains a workflow for assembling species occurrence records, integrating trait and environmental data, generating spatial predictors, fitting statistical models, producing spatial predictions, and quantifying biodiversity loss across spatial scales. The workflow is organized into sequential modules that transform raw inputs into standardized datasets, build model-ready objects, run both GLMM and Bayesian models, generate mapped outputs, and compute biodiversity metrics from posterior samples.

Intermediate outputs are saved to disk, allowing individual stages to be re-run without rebuilding the entire pipeline. Directory names reflect the main workflow components:

1_GBIF_search – Downloads and filters initial GBIF occurrence data.

2_GBIF_clean – Cleans, deduplicates, and standardizes occurrence records.

3_Geographic_range – Estimates species’ geographic ranges.

4_Traits – Compiles and formats trait data.

5_Standardized_elevation – Extracts and standardizes elevation predictors.

6_DB_dataset – Assembles the final dataset used for modelling.

7_models – Contains GLMM and Bayesian model scripts.

8_predictions – Generates spatial predictions from fitted models.

9_biodiversity_loss – Calculates biodiversity loss metrics across spatial scales.
