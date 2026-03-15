# Geographical ranges of dung beetle species/morphospecies from Colombia

This repository contains the scripts used to construct the spatial framework and species geographic ranges used in the analyses. The workflow follows the procedures described in the Range maps section of the Methods in the main manuscript and the Supplemental Methods section "*Geographical ranges of dung beetle species/morphospecies from Colombia*".

Occurrence records used to estimate species ranges were compiled from the Global Biodiversity Information Facility (GBIF) and curated taxonomic sources. Geographic ranges were approximated using the extent of occurrence (EOO) derived from occurrence records and constrained using species-specific elevational limits from a digital elevation model.

# Scripts

## Definition of the Colombian mainland boundary [1_mainland.R](./1_mainland.R) 

Defines the mainland boundary of Colombia used as the spatial mask for subsequent analyses. This step removes offshore islands and ensures that all spatial operations are performed within a consistent national boundary.

## Delineation of topographic regions [2_topographic_units.R](./2_topographic_units.R) 

Constructs the major topographic units used to structure the spatial analyses. These units are derived from HydroSHEDS drainage basins combined with mountain system boundaries, allowing the separation of major Andean slopes and surrounding lowland regions.

## Estimation of species geographic ranges [3_geographic_ranges.R](./3_geographic_ranges.R) 
Purpose

Generate geographic range maps for dung beetle species and morphospecies from Colombia based on occurrence records from [2_GBIF_clean](../2_GBIF_clean) and altitudinal constraints, following the procedures described in the main manuscript (Range maps) and the Supplemental Methods (Geographical ranges of dung beetle species/morphospecies from Colombia).

Workflow

1) Occurrence filtering: Occurrence coordinates are filtered to retain spatially unique records for each species.

2) Extent of occurrence (EOO): A geographic envelope is constructed from the occurrence points:

- 1 record: buffered point
- 2 records: buffered line connecting both points
- ≥3 records: convex hull polygon

3) Buffering of geographic envelope: The resulting geometry is buffered to account for georeferencing uncertainty and potential unsampled areas near known occurrences.

4) Restriction to study area: The buffered envelope is intersected with the Colombian mainland boundary to restrict species ranges to the study area.

5) Altitudinal constraint: Minimum and maximum elevations from species records are expanded by ±100 m to account for measurement uncertainty. A digital elevation model is used to retain only cells falling within this adjusted elevation range.

6) Range polygon generation: Elevationally suitable cells are converted to polygons representing the potential geographic range of each species.

These maps integrate both geographic and topographic constraints and provide the spatial basis for the diversity analyses presented in the manuscript.
