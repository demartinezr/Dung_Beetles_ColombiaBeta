# GBIF Occurrence Data Retrieval (GBIF_records_search.R)

Purpose

Retrieve georeferenced GBIF records for all species and subspecies in the Scarabaeinae abundance dataset, including known synonyms, to maximize coverage for constructing range maps. This script supports the "*range maps*" analyses described in the Methods and the Supplementary Online Material (SOM) “*Geographical ranges of dung beetle species/morphospecies from Colombia*”.

Workflow:

- Extract species and subspecies names from the database.
- Retrieve known synonyms via GBIF Backbone Taxonomy.
- Consolidate valid names, subspecies, and synonyms into a single list.
- Download georeferenced occurrence records from GBIF (rgbif, limit 15,000 per taxon).
- Clean data, update names to accepted taxonomy, remove duplicates.
- Produce a harmonized dataset ready for range mapping.
