# Tropical insect biodiversity loss from forest–pasture conversion is substantially underestimated across spatial scales

This repository contains the analytical workflow used to assemble biodiversity datasets, integrate environmental and trait information, fit species abundance models, generate spatial predictions, and quantify biodiversity responses to forest–pasture conversion across spatial scales in Colombia.

The pipeline combines species occurrence records, functional trait data, and environmental predictors to model dung beetle abundances using generalized linear mixed models (GLMM) and Bayesian hierarchical models. Model predictions are then used to estimate biodiversity responses under forest and pasture scenarios across spatial grains ranging from local grid cells (2 × 2 km) to ecoregional and near-national extents.

The workflow is organized into sequential modules. Each module transforms the outputs of the previous stage into standardized datasets or derived products used in subsequent analyses. Intermediate outputs are saved to disk, allowing individual stages to be re-run without rebuilding the entire pipeline.

## Workflow structure

The repository is organized into the following modules:

| Step | Module | Description |
|-----|-----|-----|
| 1 | [1_GBIF_search](./1_GBIF_search) | Download and initial filtering of GBIF occurrence records |
| 2 | [2_GBIF_clean](./2_GBIF_clean) | Cleaning, deduplication, and taxonomic standardization of occurrence records |
| 3 | [3_Geographic_range](./3_Geographic_range) | Estimation of species geographic ranges |
| 4 | [4_Traits](./4_Traits) | Compilation and formatting of functional trait data |
| 5 | [5_Standardized_elevation](./5_Standardized_elevation) | Extraction and standardization of elevation predictors |
| 6 | [6_DB_dataset](./6_DB_dataset) | Assembly of the final modelling dataset |
| 7 | [7_models](./7_models) | GLMM and Bayesian hierarchical abundance models |
| 8 | [8_predictions](./8_predictions) | Spatial predictions of species abundances under forest and pasture scenarios |
| 9 | [9_biodiversity_loss](./9_biodiversity_loss) | Quantification of biodiversity responses across spatial scales |

Each directory contains the scripts required for that stage of the workflow together with documentation describing the analyses implemented in that module.
