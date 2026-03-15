# Biodiversity loss in simulated scenearios (forest - pasture)

This directory contains scripts used to quantify biodiversity responses to forest–pasture conversion across spatial scales using species-level abundance predictions derived from the hierarchical model.

The analyses implement the sensitivity-based biodiversity metrics described in the Methods sections “*Spatial scaling, abundance-based sensitivity metrics, and species inclusion criteria*”.  “*Quantifying biodiversity loss across spatial scales using sensitivity to land-use conversion*” and "*Assessing the potential underestimation of biodiversity loss from local-scale assessments*". Predicted species abundances under forest and pasture scenarios are compared to derive forest–pasture abundance ratios, which are then summarized at multiple spatial grains: local grid cells (2 × 2 km), ecoregions, and the near-national extent.

Species inclusion in each assemblage is controlled by minimum predicted abundance thresholds (1, 3, and 10 individuals), allowing sensitivity analyses of species pool definitions.

# Scripts

## [1_db_compute_loss.R](./1_db_compute_loss.R)

Defines the core functions used to compute biodiversity metrics from predicted species abundances under forest and pasture scenarios. These functions calculate forest–pasture abundance ratios, apply species inclusion thresholds, and summarize community-level responses at different spatial aggregations (grid cells, ecoregions, and near-national scale).

Species inclusion in each assemblage is controlled by minimum predicted abundance thresholds (1, 3, and 10 individuals), allowing sensitivity analyses of species pool definitions.

The main functions implemented in this script are:

- **`avg_ratio()`** – Computes species-level forest–pasture abundance ratios and summarizes their distribution (mean, median, and quantiles of log-ratios) within a spatial unit.

- **`get_avg_cell_ratios()`** – Applies the ratio calculations at the local scale (2 × 2 km grid cells) using predicted abundances for each species.

- **`get_regional_ratios()`** – Aggregates predicted abundances across multiple grid cells to compute sensitivity metrics at broader spatial scales (ecoregions or near-national extent).

- **`get_sample_percent_decline()`** – Calculates richness-based biodiversity responses by estimating the proportion of species whose predicted abundance is higher in forest than in pasture (forest-pasture ratio > 1).

## [2_db_decline.R](./2_db_decline.R)

Implements the richness-based biodiversity loss metric used in the manuscript. Using posterior predictions from the abundance model and functions defined in db_compute_loss, this script calculates the proportion of species classified as declining (forest–pasture abundance ratio > 1) across posterior draws.

The script performs two complementary analyses across spatial scales:

1) Richness-based biodiversity loss across scales (Fig. 2a).
The first analysis compares the median proportion of declining species at three spatial scales—local (2 × 2 km grid cells), ecoregional, and near-national—to quantify how biodiversity loss changes with spatial aggregation.

2) Spatial scaling of richness-based biodiversity loss across ecoregions (Fig. 3).
The second analysis estimates the posterior density distributions of biodiversity loss for each ecoregion and compares them to the near-national estimate. These distributions illustrate how biodiversity loss varies geographically and whether regional patterns exceed local-scale losses.

## [3_db_biodiversity_loss.R](./3_db_biodiversity_loss.R)

This script implements abundance-based biodiversity loss analyses for dung beetle communities across forest–pasture landscapes. It focuses on quantifying species-level and community-level sensitivity to forest conversion and assessing how biodiversity loss varies across spatial scales and ecoregions.

Main functionalities:

1) Summarizes these metrics to evaluate how biodiversity responses change with spatial aggregation

      - Summarizes community sensitivity to forest–pasture conversion across local (2 × 2 km), ecoregional, and near-national scales.
    
      - Evaluates how median forest/pasture abundance ratios change with spatial aggregation, highlighting scale-dependent patterns of biodiversity loss (Fig. 2b).

2) variation in sensitivity among ecoregions and relative to the near-national benchmark

      - Examines variation in sensitivity among ecoregions across the species response distributions (25th, 50th, 75th percentiles).
    
      - Compares ecoregion-level responses to the near-national benchmark to identify which regions experience stronger or weaker biodiversity loss (Fig. 4a, d, g).

3) Severity of biodiversity loss at the near-national scale relative to individual ecoregions

      - Quantifies deviations of ecoregion-level sensitivity from the near-national estimate.
    
      - Captures species-level variation and uncertainty, providing a measure of how representative individual ecoregions are of broader-scale biodiversity loss (Fig. 4b, e, h).

4) Regional pooling to assess underestimation by local-scale reliance

      - Assesses potential underestimation of biodiversity loss when analyses rely on incomplete spatial coverage.
    
      - Sequentially accumulates ecoregions to evaluate cumulative sensitivity and compare pooled estimates to the near-national benchmark (Fig. 4c, f, i).
