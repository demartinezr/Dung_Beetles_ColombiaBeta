# Elevation ranges and standardized elevation

This module compiles elevation records for dung beetle species and morphospecies and derives species-level elevation ranges used in subsequent analyses.

The procedures implemented here correspond to the elevation predictor described in the Environmental predictors section of the Methods in the main manuscript, and support the analyses presented in the Supplemental Methods.

The objective of this module is to estimate species-level elevation ranges from available occurrence records and derive a standardized elevation gradient used as an environmental predictor in the Bayesian hierarchical model of species abundance.

## Workflow

Elevation values were extracted from the ALOS PALSAR Digital Elevation Model (30 m resolution). For each species or morphospecies, the following summary statistics were calculated from unique geographic records:

- minimum elevation
- maximum elevation
- mean elevation

When only a single occurrence record was available, the elevation range was expanded by ±100 m to account for limited sampling and potential measurement uncertainty.

These elevation values were then used to generate a standardized elevation variable, scaled so that −1 represents the minimum recorded elevation, 0 the mean elevation, and 1 the maximum elevation for each species. This standardized gradient is used as an environmental predictor in the Bayesian hierarchical model to estimate species-specific responses to elevation and land-use change.
