# GBIF Occurrence Data Cleaning and Integration [GBIF_records_clean.R](./GBIF_records_clean.R)

Purpose
This script cleans, standardizes, and integrates species records from multiple sources to produce a unified georeferenced dataset suitable for generating accurate species range maps. It is associated with the Methods section “*Range maps*” and the Supplementary Online Material (SOM) titled “*Geographical ranges of dung beetle species/morphospecies from Colombia*”.

Data sources:

- GBIF – Global Biodiversity Information Facility.
- IAvH-E – Colección entomológica del Instituto Humboldt.
- MEFLG – Museo Entomológico #Francisco Luis Gallego#, Universidad Nacional de Colombia.
- ICN – Instituto de Ciencias Naturales, Universidad Nacional de Colombia.
- UPTC – Museo de historia natural “Luis Gonzalo Andrade”, Universidad Pedagógica y Tecnológica de Tunja.
- CEUN-PSO – Colección entomológica de la Universidad de Nariño.
- Scarabaeinae project dataset – Field and abundance data compiled for the study.

Workflow:

1) Load species occurrence records from GBIF and entomological collections (IAvH-E, MEFLG, ICN, UPTC, CEUN-PSO).

2) Standardize species names using GBIF Backbone Taxonomy.

3) Clean data by removing duplicates, erroneous coordinates, and applying targeted corrections.

4) Incorporate validated morphospecies records.

5) Merge all cleaned records into a single georeferenced dataset ready for range mapping.
