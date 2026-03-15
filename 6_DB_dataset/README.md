# Dung beetle unified dataset

This script [db_format_for_analysis.R](./db_format_for_analysis.R) constructs the unified dataset used in the Bayesian hierarchical analysis of dung beetle abundance across Colombia.

The script compiles and harmonizes multiple data sources, integrating abundance records, sampling site metadata, species geographic ranges, biogeographic information, environmental predictors, and functional traits. The resulting dataset constitutes the final input for the statistical model described in the " Bayesian species-specific abundance modeling and prediction" section of the main  and "Full Bayesian hierarchical negative binomial mixed-effects model specification and validation" in SOM.

## Workflow

The script builds the dataset through a sequential data integration process.

First, the national dung beetle abundance database is imported and cleaned. Records from multiple sampling projects are standardized and organized into a common structure of species × sampling point × sampling day observations.

A zero-filling procedure is then applied to generate all possible species–point–day combinations, assigning zero abundance to combinations that were sampled but where a species was not detected. This ensures that non-detections are explicitly represented in the dataset.

Next, spatial metadata for sampling points are incorporated, including habitat classification, elevation values derived from the ALOS PALSAR digital elevation model, and spatial clustering information. Topographic regions and mountain slopes are assigned using spatial layers produced in the module [3_Geographic_range](../3_Geographic_range).

Species biogeographic information is then integrated to describe their distribution across Colombian regions and Andean slopes.

The script subsequently evaluates feasible species–site combinations by intersecting sampling points with species geographic range polygons derived in the module [3_Geographic_range](../3_Geographic_range). This step identifies whether a sampling location falls within the potential distribution of each species.

Environmental covariates used in the statistical model are then calculated. These include the standardized elevation gradient derived from species elevation ranges [5_Standardized_elevation](../5_Standardized_elevation) and the distance from sampling points to species range boundaries, calculated from geographic range polygons.

Functional trait information describing species ecological strategies is then merged with the dataset, including behavioral traits and morphometric traits compiled in module [4_traits](../4_Traits).

Finally, abundance data are aggregated across sampling days, habitats not included in the analysis are removed, and species–site combinations outside geographic ranges are filtered. The resulting dataset represents the complete set of observations used in the Bayesian hierarchical model of dung beetle abundance.

Final dataset

The final dataset links dung beetle abundance with environmental predictors, geographic constraints, and functional trait information. This dataset constitutes the input used to estimate species responses to elevation, land-use change, and functional traits in the Bayesian hierarchical model presented in the manuscript.
