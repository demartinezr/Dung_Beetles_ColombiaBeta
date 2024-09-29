# Directory
    setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata/GBIF_data")
# Packages
    library(readxl)
    library(rgbif)
    library(dplyr)
    library(ggplot2)
    library(rnaturalearth)
    library(rnaturalearthdata)
    library(sf)

# Dung beetles Dataset #
    IBDdata <- read_excel("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/Scarabaeinae_database_2024.xlsx", sheet="Scarabaeinae_database_2024")
#
# List of species and subspecies #
    species_list <- data.frame(subset(IBDdata, taxonRank == "Species"))
    species_list <- sort(unique(species_list$scientificName)) 
    species_list <- gsub("_", " ", species_list)
    subspecies_list <- data.frame(subset(IBDdata, taxonRank == "Subspecies"))
    subspecies_list <- sort(unique(subspecies_list$scientificName)) 
    subspecies_list <- gsub("_", " ", subspecies_list)
    species_list <- sort(c(species_list, subspecies_list))
# List of synonyms
    synonym_search <- list()
    for (sp in species_list) {
      result<- name_lookup(query = sp, rank = c("species", "subspecies"), status = "SYNONYM")$data
      if (!is.null(result)) {
        result$scientificName1 <- sp
        synonym_search[[sp]] <- result
      }
    }
    synonym_search <- lapply(synonym_search, function(x) {
      all_columns <- unique(unlist(lapply(synonym_search, colnames)))
      missing_columns <- setdiff(all_columns, colnames(x))
      x[missing_columns] <- NA  
      return(x)
    })
    all_synonyms <- do.call(rbind, synonym_search)
    rownames(all_synonyms) <- NULL
    #saveRDS(all_synonyms, "all_synonyms")
    #write.table(all_synonyms, "all_syonyms.txt")
    synonym_list <- unique(sort(all_synonyms$canonicalName))
# Species, morfospecies and synonyms list
    all_taxa <- sort(unique(c(species_list, synonym_list)))
#
# GBIF search #
    result <- occ_data(scientificName = all_taxa, limit=15000, hasCoordinate = TRUE)
    data_list <- lapply(result, function(x) x$data)
    records <- data.frame(bind_rows(data_list, .id = "scientificName"))
    list_cols <- sapply(records, is.list)  #convert columns list in characters
    records[list_cols] <- lapply(records[list_cols], function(x) sapply(x, toString))
    #saveRDS(records, "GBIF_2024-09-21.rds")
    #write.table(records, "GBIF_2024-09-21.txt", row.names = FALSE, sep=",")
    #
  # Update species name by GBIF Backbone Taxonomy
    canonical_names <- all_synonyms$canonicalName[match(records$scientificName, all_synonyms$canonicalName)]
    accepted_names <- all_synonyms$scientificName1[match(records$scientificName, all_synonyms$canonicalName)]
        # Actualizar la columna scientificName1 en records
    records <- records %>%
      mutate(
        scientificName1 = if_else(!is.na(canonical_names), accepted_names, scientificName)
      )
    records$coordinates <- paste(records$decimalLatitude, records$decimalLongitude, sep="_")
    records <- records[!is.na(records$coordinates) & !duplicated(records[c("species", "coordinates")]), ]
    #
# Check records #
#
    species_count <- as.data.frame(table(sort(records$species)))
    #names(species_count) <- c("species", "count")
# subset records of Deltochilum carinatum
    sp <- data.frame(subset(records, scientificName1 == "Ateuchus irinus"))
#
# base maps #
    south_america <- ne_countries(scale = "medium", continent = "South America", returnclass = "sf")
    central_america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
    central_america <- central_america[central_america$name %in% c("Belize", "Costa Rica", "El Salvador", "Mexico", "Guatemala", "Honduras", "Nicaragua", "Panama"), ]
    ecoreg <- st_read(dsn = "C:/Users/Dell-PC/Dropbox/CO_DBdata/SIG/Ecoregions2017", layer = "eco_col") #Pais
    ecoreg <- st_transform(ecoreg, st_crs("+proj=longlat +datum=WGS84"))
    col <- ne_states(country = "colombia", returnclass = "sf")
# plot 
    ggplot() +
     geom_sf(data = south_america, fill = "lightgrey", color = "black") +
     geom_sf(data = central_america, fill = "lightgrey", color = "black") +
     geom_point(data = sp, aes(x = decimalLongitude, y = decimalLatitude), color = "blue")
