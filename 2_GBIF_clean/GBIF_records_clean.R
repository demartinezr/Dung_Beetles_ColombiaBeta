# This script cleans, standardizes, and integrates species records from multiple 
# sources, including GBIF and biological collections from IAvH-E, MEFLG, ICN, UPTC, 
# and CEUN-PSO, to produce a unified georeferenced dataset. It updates species names 
# based on the GBIF Backbone Taxonomy, removes duplicated and erroneous coordinates, 
# incorporates validated morphospecies records, and applies targeted corrections for 
# specific taxa. The script then merges cleaned species and morphospecies occurrences 
# from all sources into a final consolidated dataset suitable for generating accurate 
# species range maps.

setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata/GBIF_data")
# packages
  library(dplyr)
  library(rgdal)
  library(sf)
  library(readxl)
#
# Species data: Records from GBIF and entomological collections
  #  
  # GBIF dataset
    all_synonyms <- readRDS("all_synonyms.rds")    
    GBIF <- readRDS("GBIF_2024-09-21.rds")
    # Updated species names by GBIF Backbone Taxonomy
    canonical_names <- all_synonyms$canonicalName[match(GBIF$scientificName, all_synonyms$canonicalName)]
    accepted_names <- all_synonyms$scientificName1[match(GBIF$scientificName, all_synonyms$canonicalName)]
    GBIF <- GBIF %>%
      mutate(
      scientificName1 = if_else(!is.na(canonical_names), accepted_names, scientificName)
      )
    GBIF$coordinates <- paste(GBIF$decimalLongitude, GBIF$decimalLatitude, sep="_")
    GBIF <- GBIF[!is.na(GBIF$coordinates) & !duplicated(GBIF[c("scientificName1", "coordinates")]), ]
    GBIF <- GBIF[!(GBIF$scientificName1=="Ontherus sulcator"), ]
    #
  # Entomological collections dataset (IAvH-E, MEFLG, UPTC, CEUN-PSO, Master degree 
  # DEMR data from Orinoquía IAvH-E, ICN, CECC-CALT)
    SIB_COL <- read_excel("DB-Col_2024-01-23.xlsx", sheet = "Conjunto-2020-07-23")
      SIB_COL$coordinates <- paste(SIB_COL$decimalLongitude, SIB_COL$decimalLatitude, sep = "_")
      SIB_COL$scientificName1 <- paste(SIB_COL$identificationRemarks)
      SIB_COL <- SIB_COL[!is.na(SIB_COL$coordinates) & !duplicated(SIB_COL[c("scientificName1", "coordinates")]), ]
  #
  # Dung beetles data set project
    IBD_data <- read_excel("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/Scarabaeinae_database_2024.xlsx", sheet = "Scarabaeinae_database_2024")
      IBD_data$scientificName1 <- gsub("_", " ", IBD_data$scientificName)
      IBD_data$decimalLongitude <- as.numeric(IBD_data$lon_all_points)
      IBD_data$decimalLatitude <- as.numeric(IBD_data$lat_all_points)
      IBD_data$coordinates <- paste(IBD_data$lon_all_points, IBD_data$lat_all_points, sep = "_")
      IBD_data <- IBD_data[!(IBD_data$abundance =="NA") & !duplicated(IBD_data[c("scientificName", "coordinates")]), ]
      IBD_data <- IBD_data[!(IBD_data$abundance == 0), ]
    #
    # cleaning
  records <- GBIF
  records <- records[records$coordinates != "-74.266667_4.05", ] # Delete records of A. villosus, C. juvencus, C. telamon and P. haroldi from páramo
    at_cereus <- data.frame(subset(SIB_COL, scientificName1=="Ateuchus sp. 08H")) # Select records of A. cereus from entomological collections data base
    at_cereus <- mutate(at_cereus, scientificName1 = ifelse(scientificName1 == "Ateuchus sp. 08H", "Ateuchus cereus", scientificName1))
  records <- merge(records, at_cereus, by = intersect(names(records), names (at_cereus)), all = TRUE) # add records of A. cereus from Vichada
    at_irinus <- data.frame(subset(SIB_COL, identificationRemarks == "Ateuchus irinus")) # Select records of A. irinus of Nariño from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, at_irinus, by = intersect(names(records), names(at_irinus)), all = TRUE) # add records for A. irinus
  records <- records[!(records$scientificName1=="Ateuchus irinus" & records$coordinates %in% c("-5.983333_-2.416667", "-74.855_2.7975")),] # Delete A. irinus records from Ocean
    at_murrayi <- data.frame(subset(SIB_COL, scientificName1=="Ateuchus sp. 09H")) # select records for A. murrayi from Collections data base
    at_murrayi <- mutate(at_murrayi, scientificName1= ifelse(scientificName1=="Ateuchus sp. 09H", "Ateuchus murrayi", scientificName1)) # Updated morphospiecies to species
  records <-merge(records, at_murrayi, by = intersect(names(records), names (at_murrayi)), all = TRUE) # Add records for A. murrayi from Vichada
  records <- records[!(records$scientificName1=="Ateuchus murrayi" & records$coordinates == "-5.983333_-2.416667"), ]  
    at_substriatus <- data.frame(subset(SIB_COL, identificationRemarks == "Ateuchus sp. 10H")) # Select records of Ateuchus substriatus from entomological collections data base
    at_substriatus <- mutate(at_substriatus, scientificName1 = ifelse(scientificName1=="Ateuchus sp. 10H", "Ateuchus substriatus", scientificName1)) #  Updated morphospiecies to species
  records <- merge(records, at_substriatus, by = intersect(names(records), names(at_substriatus)), all = TRUE) # add records for A. substriatus from Vichada
    cdm_gers <- data.frame(subset(IBD_data, scientificName1 =="Canthidium gerstaeckeri"))
  records <-  merge(records, cdm_gers, by=intersect(names(records), names(cdm_gers)), all = TRUE)
    cdm_on <- data.frame(subset(SIB_COL, identificationRemarks == "Canthidium sp. 36H")) # Select records of Canthidium onitoides from "morphospecies unification SIB "DB-Col-2022-07-23"
    cdm_on <- mutate(cdm_on, scientificName1 = ifelse(scientificName1 == "Canthidium sp. 36H", "Canthidium onitoides", scientificName))  # replace names
  records <- merge(records, cdm_on, by = intersect(names(records), names(cdm_on)), all = TRUE) # add records for C onitoides
    ctn_fu <- subset(SIB_COL, identificationRemarks =="Canthon fulgidus")# Select records of Canthon fulgidus from "morphospecies unification with SIB "DB-Col-2022-07-23"
    ctn_fu <- mutate(ctn_fu, scientificName1 = ifelse(scientificName1 == "Canthon fulgidus", "Canthon fulgidus martinezi", scientificName1)) # replace names
  records <- merge(records, ctn_fu, by = intersect(names(records), names(ctn_fu)), all=TRUE)# add records for C. fulgidus martinezi
    c_juv <- data.frame(subset(SIB_COL, identificationRemarks == "Canthon juvencus")) # Select records of Canthon juvencus from "morphospecies unification SIB "DB-Col-2022-07-23" our data
  records <- merge(records, c_juv, by = intersect(names(records), names(c_juv)), all=TRUE)# add records for C. juvencus
  records <- records[!(records$scientificName1=="Canthon juvencus" & records$coordinates %in% 
                    c("-75078_10823", "-74907_9027", "-74831_4936", "-74.84_9741",
                      "-75.07859583_6.913268902", "-74.25_4.05", "NA_NA")),] # Remove atypical records of C. juvencus
    ctn_li <- data.frame(subset(SIB_COL, scientificName1 == "Canthon lituratus")) # Select records of Canthon lituratus from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, ctn_li, by = intersect(names(records), names(ctn_li)), all = TRUE) # add records for C lituratus
  records <- records[!(records$scientificName1=="Canthon lituratus" & records$coordinates %in% c("-75078_10823","-74.84_9741", "-59.83333333_-2.416666667", "-74.083333_11.116667", "-59.833333_-2.416667")),] # Remove atypical altitudinal records of C. lituratus
  records <- records[!(records$scientificName1=="Canthon luteicollis" & records$stateProvince %in% c("Arauca", "Casanare", "Vichada")),] # Remove records of C. luteicollis from Orinoquia
  records <- records[!(records$scientificName1=="Canthon luteicollis" & records$coordinates == "-77.6_-0.46667"),] # Remove atypical altitudinal records of de C. luteicollis
  records <- records[!(records$scientificName1=="Canthon pallidus" & records$coordinates =="-74.166667_2.666667"),]
  records <- records[!(records$scientificName1=="Canthon semiopacus" & records$coordinates =="-79.24529_-3.99067"),]
  records <- records[!(records$scientificName1 == "Canthon septemmaculatus" & records$decimalLatitude <= 5),]  # Remove records with latitude <= 5 for C. septemmaculatus
  records <- records[!(records$scientificName1 == "Canthon septemmaculatus" & records$coordinates %in% c("1.59858_49.06679")),]  # Remove outliers coordinates for C. septemmaculatus
  records <- records[!(records$scientificName1 == "Canthon septemmaculatus" & records$decimalLongitude >= -70),]  # Remove atypical coordinates in longitude for C. septemmaculatus, maybe it is C. septemmaculatus linearis
    ctn_su <- data.frame(subset(SIB_COL, identificationRemarks == "Canthon subhyalinus")) # Select records of Canthon subhyalinus from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, ctn_su, by = intersect(names(records), names(ctn_su)), all = TRUE) # add records for Canthon subhyalinus
  records <- records[records$coordinates != "-74831_4936", ] # Remove atypical records of C. subhyalinus
  records <- records[!(records$scientificName1=="Canthon subhyalinus" & records$coordinates =="-75.07859583_6.913268902"),] # Remove atypical altitudinal records of C. subhyalinus
  records <- records[!(records$scientificName1=="Canthon subhyalinus" & records$decimalLatitude <= 1.7),] # Remove records less than 1.7 Longitude of C. subhyalinus
  records <- records[!(records$scientificName1=="Canthon subhyalinus" & records$stateProvince == "Meta"),] # Remove records from llanos of C. subhyalinus
  records <- records[!(records$scientificName1=="Coprophanaeus corythus" & records$stateProvince %in% c("Meta", "Casanare", "Vichada")),] # remove records of C. corythus from Orinoquia
  records <- records[!(records$scientificName1=="Coprophanaeus corythus" & records$decimalLatitude >= 23),] # delete records less than 23 Latitude of C. corythus
    c_edmondsi <- data.frame(subset(SIB_COL, identificationRemarks == "Coprophanaeus edmondsi")) # Select records of C. edmondsi from Entomological collections
  records <- merge(records, c_edmondsi, by = intersect(names(records), names(c_edmondsi)), all=TRUE)# add records for C. edmondsi in GBIF C. conocephalus
  records <- records[!(records$scientificName1 == "Coprophanaeus jasius" & records$stateProvince %in% c("Arauca", "Antioquia", "Canindeyu", "Sucre", "Atlántico")), ] # remove records of C. jasius from Atlantico department
  records <- records[!(records$scientificName1=="Coprophanaeus morenoi" & records$coordinates == "-78.42578_-1.39315"),] # Remove atypical altitudinal records of C. morenoi
  records <- records[!(records$scientificName1 == "Coprophanaeus telamon" & records$stateProvince %in% c("Antioquia","Caldas", "Cesar", "Córdoba", "Cundinamarca","La Guajira", "Santander", "Sucre", "Valle del Cauca")), ] # remove some departments from the distribution of C. telamon
  records <- records[!(records$scientificName1 == "Coprophanaeus telamon" & records$county == "Puerto Boyacá"), ] # Remove C. telamon from Magdalena Valley (It really is C. corythys)
  records <- records[!(records$scientificName1 == "Coprophanaeus telamon" & records$decimalLongitude <= -80), ] # Remove C. telamon from Central America (It really is C. corythys)
    c_campbe <- data.frame(subset(SIB_COL, scientificName1 == "Cryptocanthon campbellorum")) # select records of C. campbellorum from entomological collections (Master Thesis DEMR)
  records <- merge(records, c_campbe, by=intersect(names(records), names(c_campbe)), all = TRUE) # add records of C. campbellorum
    c_mail <- data.frame(subset(SIB_COL, scientificName1 == "Cryptocanthon mailinae")) # select records of C. mailinae from entomological collections (IAvH-E)
  records <-  merge(records, c_mail, by = intersect(names(records), names(c_mail)), all = TRUE) # add records of C. mailinae
  records <- records[!(records$scientificName1 == "Deltochilum carinatum" & records$country == "Costa Rica"), ] # Remove D. carinatum from Costa Rica
  records <- records[!(records$scientificName1 == "Deltochilum carinatum" & records$coordinates %in% c("-76.1_1.333333", "-76.544889_1.243611")), ] # Remove D. carinatum from Andes
  records <- records[!(records$scientificName1=="Deltochilum guildingii" & records$stateProvince == "Caquetá"),] # remove D. guildingii from Caquetá
  records <- records[!(records$scientificName1=="Deltochilum guildingii" & records$coordinates=="-74.083333_11.116667"), ] # remove atypical altitudinal records of D. guildingii
  records <- records[!(records$scientificName1=="Deltochilum hypponum" & records$coordinates %in% c("-75.428647_0.838333", "-75.246667_2.935278", "-76.633333_1.133333", "-73.388056_4.5925", "-73.408611_5.75")), ] # remove atypical altitudinal records of D. hypponum
  records <- records[!(records$scientificName1=="Deltochilum mexicanum" & records$coordinates %in% c("-71.8375_6.513333333", "-73.38805556_4.5925", "-78.25_1.283333333", "-74.72043_2.71329")), ] # remove atypical altitudinal records of D. burmeisteri
  records <- records[!(records$scientificName1=="Deltochilum mexicanum" & records$coordinates %in% c("Guatemala","Mexico")), ] # Remove D. burmeisteri from Guatemala and Mexico 
  records <- records[!(records$scientificName1=="Deltochilum molanoi" & records$coordinates %in% c("-76.633333_2.283333", "-75.627222_5.661111")),] # remove points greater than 76 Longitude D, molanoi
  records <- records[!(records$scientificName1=="Deltochilum orbiculare" & records$decimalLongitude <= -76.4),] # remove records greater than 76 Longitude D. orbiculare
  records <- records[!(records$scientificName1=="Dichotomius agenor" & records$decimalLongitude <= -81),] # remove records greater than  81 Long D agenor
  records <- records[!(records$scientificName1=="Dichotomius agenor" & records$stateProvince %in% c("Boyacá",  "Meta","Casanare", "Vichada", "BoyacÃ¡")), ] # Remove D. agenor departments based on Montoya et al. 2021
  records <- records[!(records$scientificName1=="Dichotomius agenor" & records$coordinates %in% c("-74.84_9741", "-74907_9027", "-74831_4936", "NA_NA", "-74907_9027","-75078_10823", "-74831_4936", "-74.06785_11.0973")), ] # remove atypical records of D. agenor
    d_andresi <- data.frame(subset(SIB_COL, scientificName1=="Dichotomius andresi")) # select records of D. andresi from Entomological collections data
  records <- merge(records, d_andresi, by = intersect(names(records), names(d_andresi)), all = TRUE) # add records of D. andresi from MEFLG 
  records <- records[!(records$scientificName1=="Dichotomius andresi" & !duplicated(records$coordinates)),]
  records <- records[!(records$scientificName1=="Dichotomius belus" & records$coordinates %in% c("-49.866665_5.091945", "-71.14925_6.563528")), ] # Remove outlier records outside the country and Arauca for Dichotomius belus
  records <- records[!(records$scientificName1=="Dichotomius belus" & records$coordinates %in% c("-74.53146_6.27808", "-74.73175_5.73692", "-74.554144_6.460248", "-74.572882_6.494321")), ] # remove atypical records of  D. belus
  records <- records[!(records$scientificName1 == "Dichotomius belus" & records$stateProvince %in% c("Cesar", "Atlántico", "Magdalena")), ] # Remove D. belus departments
    d_compre <- data.frame(subset(SIB_COL, scientificName1=="Dichotomius compressicollis")) # Select records of Dichotomius compressicollis from Entomological collections dataset
  records <- merge(records, d_compre, by=intersect(names(records), names(d_compre)), all=TRUE) # add recors of D. compressicollis
    di_gan <- data.frame(subset(SIB_COL, identificationRemarks == "Dichotomius sp. 10H")) # Select records of Dichotomius gandinii from "morphospecies unification with SIB "DB-Col-2022-07-23"
    di_gan <- mutate(di_gan, scientificName1 = ifelse(scientificName1 == "Dichotomius sp. 10H", "Dichotomius gandinii", scientificName1)) # replace names
  records <- merge(records, di_gan, by = intersect(names(records), names(di_gan)), all = TRUE) # add records of D. gandinii
    di_glo <- data.frame(subset(SIB_COL, identificationRemarks == "Dichotomius globulus")) # Select records of Dichotomius globulus from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, di_glo, by = intersect(names(records), names(di_glo)), all = TRUE) # add records of D. globulus
    di_mam <- data.frame(subset(SIB_COL, identificationRemarks == "Dichotomius mamillatus")) # Select records of Dichotomius mamillatus from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, di_mam, by = intersect(names(records), names(di_mam)), all = TRUE) # add records of D mamillatus
  records <- records[!(records$scientificName1=="Dichotomius nisus" & records$decimalLatitude <= 0),] # Remove D nisus from Brazil, because we don't like include Amazonia. D nisus inhabit open areas such as Orinoquia, el Cerrado, Catinga ...
  records <- records[!(records$scientificName1=="Dichotomius nisus" & records$coordinates %in% c("-73.885261_11.12191", "-74.85237_5.01152", "-76.79082_1.07046")),] # remove outliers records of D. nisus from Caribe, Magdalena and Putumayo.
    d_ohausi <- data.frame(subset(SIB_COL, scientificName1=="Dichotomius ohausi")) # select records of D. ohausi from literature Arias & Vaz-de-Mello 2024
  records <- merge(records, d_ohausi, by=intersect(names(records), names(d_ohausi)), all=TRUE) # add records of D. ohausi
  records <- records[!(records$scientificName1 =="Dichotomius quadrilobatus" & records$coordinates =="-77_1.15"),]
  records <- records[!(records$scientificName1 == "Dichotomius quinquelobatus" & records$stateProvince %in% c("Cauca", "Valle del Cauca", "Quindío","Risaralda", "Caldas")), ] # remove D. quinquelobatus departments
  records <- records[!(records$scientificName1=="Dichotomius quinquelobatus" & records$coordinates %in% c("-78.433692_-0.367544", "-79.24529_-3.99067", "-77.21667_0.05")),] # remove atypical points of D. quinquelobatus
    d_recli <- data.frame(subset(SIB_COL, scientificName1=="Dichotomius reclinatus")) # select records of D. reclinatus from literature Arias & Vaz de Mello 2023
  records <- merge(records, d_recli, by=intersect(names(records), names(d_recli)), all=TRUE) # add records of D. reclinatus
  records <- records[!(records$scientificName1=="Dichotomius tristis" & records$coordinates == "-76.050835_5.707222"),] # remove outlier records of Dichotomius tristis from Antioquia
    di_val <- data.frame(subset(SIB_COL, identificationRemarks == "Dichotomius validipilosus")) # Select records of Dichotomius validipilosus from "morphospecies unification with SIB "DB-Col-2022-07-23" and our data
  records <- merge(records, di_val, by = intersect(names(records), names(di_val)), all = TRUE) # add records of D validipilosus
  records <- records[!(records$scientificName1=="Dichotomius validipilosus" & records$decimalLongitude >= -60),] # Remove records less than -60 Long D. validipilosus
  records <- records[!(records$scientificName1=="Digitonthophagus gazella" & records$country != "Colombia"), ] # records from Colombia
  records <- records[!(records$scientificName1=="Digitonthophagus gazella" & records$coordinates %in% c("-72.536611_11.145056", "-76.068569_3.105836")),] # remove atypical records of D. gazella
  records <- records[!(records$scientificName1=="Eurysternus caribaeus" & records$coordinates %in% c("-52.95_5.36667", "-75.7166666666667_4.11666666666667", "-72.084167_7.959444")),] # remove atypical altitudinal records of E. caribaeus
  records <- records[!(records$scientificName1=="Eurysternus cayennensis" & records$coordinates %in% c("-88.7_17.133333", "-52.95_5.36667", "-76.6833_3.00333")), ] # remove atypical altitudinal records of E. cayennensis
  records <- records[!(records$scientificName1=="Eurysternus contractus" & records$coordinates %in% c("-77.3333_6.00667", "-70.533333_0.233333", "-78.936729_-1.07175", "-79.005517_-2.900755", "-75.666553_1.508061", "-77.216667_0.05")), ] # remove atypical records of E. contractus
  records <- records[!(records$scientificName1=="Eurysternus foedus" & records$decimalLatitude <= -20),] # remove records less than -20 Lat E. foedus
  records <- records[!(records$scientificName1=="Eurysternus hamaticollis" & records$decimalLongitude <=-80),] # remove records of E. hamaticollis from centro america
  records <- records[!(records$scientificName1=="Eurysternus hamaticollis" & records$coordinates %in% c("-76.6833_3.00333")), ] # remove atypical altitudinal records of E. hamaticollis from Cauca valley
  records <- records[!(records$scientificName1=="Eurysternus hypocrita" & records$coordinates %in% c("-52.95_5.36667", "-76.6833_3.00333", "-80.7531_-1.75806")), ] # remove atypical altitudinal records of E. hypocrita
  records <- records[!(records$scientificName1 == "Eurysternus hypocrita" & records$stateProvince %in% c("Antioquia","Santander")), ] # Remove E. hypocrita from Magdalena Valley
  records <- records[!(records$scientificName1=="Eurysternus marmoreus" & records$coordinates == "-74.033333_11.333333"), ] # remove atypical altitudinal points of E. marmoreus
  records <- records[!(records$scientificName1=="Eurysternus plebejus" & records$coordinates %in% c("19.083333_-32.079167","19.253889_-31.812778","-79.195752_-4.004985", "-79.199326_-3.986796")), ] # remove atypical altitudinal records of E. plebejus
  records <- records[!(records$scientificName1=="Eurysternus squamosus" & records$coordinates == "-76.1_1.333333"),] # remove outliers altitudinal records of E. squamosus
  records <- records[!(records$scientificName1=="Eurysternus wittmerorum" & records$coordinates %in% c("-79.199326_-3.986796", "-76.6833_3.00333")),]
  records <- records[!(records$scientificName1=="Eurysternus wittmerorum" & records$decimalLatitude >= 5),] # Remove E. wittmerorum from Boyaca
    G_aeru <- data.frame(subset(SIB_COL, identificationRemarks == "Gromphas aeruginosa")) # import records by literature for G. aeruginosa Cupello & Vaz de Mello 2013
  records <- merge(records, G_aeru, by = intersect(names(records), names(G_aeru)), all = TRUE) # add records of D validipilosus
  records <- records[!(records$scientificName1 == "Gromphas aeruginosa" & records$decimalLongitude >= -40),] # Remove G. aeruginosa outliers
  records <- records[!(records$scientificName1 == "Gromphas aeruginosa" & records$decimalLatitude >= 4),] # Remove G. aeruginosa outliers
  records <- records[!(records$scientificName1=="Gromphas lemoinei" & records$coordinates == "-74.081944_4.61"),] # Remove records of G lemoinei from Bogota DC ...
    H_ach <- data.frame(subset(SIB_COL, identificationRemarks == "Homocopris achamas")) # Select records from Nariño by "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, H_ach, by = intersect(names(records), names(H_ach)), all =T) # add records of H. achamas
  records <- records[!(records$scientificName1=="Homocopris achamas" & duplicated(records$coordinates)), ] # remove duplicated records
  records <- records[!(records$scientificName1=="Homocopris achamas" & records$stateProvince %in% c("Quindío", "Risaralda")), ] # remove duplicated records
    o_brev <- data.frame(subset(SIB_COL, identificationRemarks == "Ontherus brevicollis")) # Select records of Ontherus brevicollis from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, o_brev, by = intersect(names(records), names(o_brev)), all = TRUE) # add records of O. brevicollis from Nariño by "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- records[!(records$scientificName1=="Ontherus brevicollis" & duplicated(records$coordinates)), ] # remove duplicated records
  records <- records[!(records$scientificName1=="Ontherus brevicollis" & records$coordinates == "-73.38805556_4.5925"), ] # remove outliers altitudinal records for O. brevicollis
  records <- records[!(records$scientificName1=="Ontherus diabolicus" & records$coordinates == "-72.683333_5.433333"),] # remove atypical altitudinal records of O. diabolicus
  records <- records[!(records$scientificName1=="Ontherus diabolicus" & records$stateProvince=="Huila"),] # remove outlier records of O. diabolicus from Huila  
    on_inc <- data.frame(subset(SIB_COL, identificationRemarks == "Ontherus incisus")) # Select records of Ontherus incisus from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, on_inc, by = intersect(names(records), names(on_inc)), all = TRUE) # add records of O. incisus from Nariño
  records <- records[!(records$scientificName1=="Ontherus kirschii" & records$coordinates %in% c("-72.840833_7.438889", 	"-73.388056_4.5925")),] # remove atypical altitudinal records of O. kirschii
    on_lun <- data.frame(subset(SIB_COL, identificationRemarks == "Ontherus lunicollis")) # Select records of Ontherus lunicollis from Nariño by "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, on_lun, by = intersect(names(records), names(on_lun)), all = TRUE) #a add records of O. lunicollis
  records <- records[!(records$scientificName1=="Ontherus lunicollis" & duplicated(records$coordinates)), ] # remove duplicated records
  records <- records[!(records$scientificName1=="Ontherus lunicollis" & records$coordinates %in% c("-75.246667_2.935278", "-72.61869444_6.66511111", "-75.495361_4.469028")), ] # remove atypical altitudinal records of O. lunicollis
  records <- records[!(records$scientificName1=="Ontherus pubens" & records$coordinates == "-78.42578_-1.39315"),] # remove atypical altitudinal records of O. pubens
  records <- records[!(records$scientificName1=="Onthophagus acuminatus" & records$decimalLongitude >= -72.7),] # remove records less than -73 Lat O acuminatus
  records <- records[!(records$scientificName1=="Onthophagus acuminatus" & records$coordinates %in% c("-73.1_9.85", "-73.458333_5.746667", "-76.068569_3.105836", "NA_NA")),] # remove outliers altitudinal records of O. acuminatus
    o_bide <- data.frame(subset(SIB_COL, identificationRemarks == "Onthophagus bidentatus")) # Select records of Onthophagus bidentatus from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, o_bide, by = intersect(names(records), names(o_bide)), all = TRUE) # add records of O. bidentatus from IAvH collection, localities Nariño by SIB "DB-Col-2022-07-23"
  records <- records[!(records$scientificName1=="Onthophagus bidentatus" & duplicated(records$coordinates)), ]
  records <- records[!(records$scientificName1=="Onthophagus bidentatus" & records$coordinates %in% c("-74907_9027", "-74.84_9741", "-75078_10823", "NA_NA", "-14.5_14.5")), ] # remove atypical altitudinal records of O. bidendatus
    o_curv <- data.frame(subset(SIB_COL, identificationRemarks == "Onthophagus curvicornis"))
  records <- merge(records, o_curv, by = intersect(names(records), names(o_curv)), all = TRUE) # add records of O. curvicornis from Nariño by "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- records[!(records$scientificName1 == "Onthophagus curvicornis" & duplicated(records$coordinates)),]
  records <- records[!(records$scientificName1=="Onthophagus curvicornis" & records$decimalLatitude >= 10),] # remove records greater than 10 Lat O curvicornis
  records <- records[!(records$scientificName1=="Onthophagus curvicornis" & records$decimalLongitude >= -70),] # remove records gerater than 70 Long O curvicornis
  records <- records[!(records$scientificName1=="Onthophagus curvicornis" & records$coordinates %in% c("-72_4", "-71.86113_6.47894", "-72.309928_5.608447", "-72.987472_5.946583")),] # remove atypical altitudinal records of O. curvicornis
  records <- records[!(records$scientificName1 == "Onthophagus marginicollis" & records$decimalLongitude >= -50), ] # remove atypical records of O. marginicollis
    o_tran <- data.frame(subset(SIB_COL, identificationRemarks == "Onthophagus transisthmius")) # Select records of Onthophagus transisthmius from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, o_tran, by = intersect(names(records), names(o_tran)), all = TRUE) # Add records of O. transithmius
    ox_ebe <- data.frame(subset(SIB_COL, identificationRemarks == "Oxysternon ebeninum")) # Select records of Oxysternon ebeninum from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, ox_ebe, by = intersect(names(records), names(ox_ebe)), all = TRUE) # add records of O. ebeninum, our data Putumayo
    ox_sil <- data.frame(subset(SIB_COL, identificationRemarks == "Oxysternon silenus")) # Select records of Oxysternon silenus from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, ox_sil, by = intersect(names(records), names(ox_sil)), all = TRUE) # add records of O. silenus
  records <- records[!(records$scientificName1=="Phanaeus bispinus" & records$coordinates=="-79.199939_-3.982992"),] # remove outlier records from Western Cordillera Ecuador
  records <- records[!(records$scientificName1 == "Phanaeus chalcomelas" & records$coordinates %in% c("-78.433692_-0.367544", "-79.24529_-3.99067", "-80.74694_-1.78167")), ] # remove atypical altitudinal records P. chalcomelas
  records <- records[!(records$scientificName1=="Phanaeus haroldi" & records$decimalLongitude >= -65),] # remove records less than -60 Long of P. haroldi
  records <- records[!(records$scientificName1=="Phanaeus haroldi" & records$coordinates %in% c("-72.691667_5.434722", "-72.691667_5.434722")),] # remove atypical altitudinal records of P. haroldi
  records <- records[!(records$scientificName1=="Phanaeus meleagris" & records$coordinates =="-79.199326_-3.986796"),]
  records <- records[!(records$scientificName1 == "Phanaeus pyrois" & records$country != "Colombia"), ] #extraer solo información de P. pyrois para Colombia
  records <- records[!(records$scientificName1 == "Pseudocanthon perplexus"),] # remove records from P. perplexus, and we work with unification with SIB "DB-Col-2022-07-23" 
    ps_per <- data.frame(subset(SIB_COL, identificationRemarks == "Pseudocanthon sp. 01H")) # Select records of Pseudocanthon perplexus (sp 01H) from "morphospecies unification with SIB "DB-Col-2022-07-23"
    ps_per <- mutate(ps_per, scientificName1 = ifelse(scientificName1 == "Pseudocanthon sp. 01H", "Pseudocanthon perplexus", scientificName1)) # replace names
  records <- merge(records, ps_per, by = intersect(names(records), names(ps_per)), all = TRUE)# add records of P. perplexus
  records <- records[!(records$scientificName1== "Pseudocanthon perplexus" & records$decimalLatitude== "9027"),] # remove atypical altitudinal records of Pseudocanthon perplexus
  records <- records[!(records$scientificName1 == "Pseudocanthon xanthurus"),] # remove records from P. xanthurus, and we work with unification with SIB "DB-Col-2022-07-23" 
    ps_xan <- data.frame(subset(SIB_COL, identificationRemarks == "Pseudocanthon sp. 02H")) # Select records of Pseudocanthon xanthurus (sp 02H) from "morphospecies unification with SIB "DB-Col-2022-07-23"
    ps_xan <- mutate(ps_xan, scientificName1 = ifelse(scientificName1 == "Pseudocanthon sp. 02H", "Pseudocanthon xanthurus", scientificName1)) # replace names
  records <- merge(records, ps_xan, by = intersect(names(records), names(ps_xan)), all = TRUE)# add records of P. xanthurus
  records <- records[!(records$scientificName1=="Scatimus ovatus" & records$decimalLongitude >= -72.4),] # remove records less than -72.4 Long S. ovatus
  records <- records[!(records$scientificName1=="Scatimus ovatus" & records$decimalLatitude <= 4.9),] # remove records less than 4.9 Lat S. ovatus
  records <- records[!(records$scientificName1=="Scybalocanthon kelleri" & records$coordinates=="-73.595_3.045556"), ]# remove atypical altitudinal records of S. kelleri
    sc_dar <- data.frame(subset(SIB_COL, scientificName == "Scybalocanthon darlingtoni"))# Select records of Scybalocanton darlingtoni from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, sc_dar, by = intersect(names(records), names(sc_dar)), all = TRUE)# add records of S. darlingtoni by Silva & Valois 2019
  records <- records[!(records$scientificName1=="Scybalocanthon darlingtoni" & records$coordinates == "-74.166664_10.916667"),] # remove records greater than 11 Lat S. darlingtoni
    sc_mar <- data.frame(subset(SIB_COL, scientificName == "Scybalocanthon martinezi"))# Select records of S. martinezi from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, sc_mar, by = intersect(names(records), names(sc_mar)), all = TRUE)# add recrods of S. martinezi by Silva & Valois 2019
    sc_sex <- data.frame(subset(SIB_COL, identificationRemarks == "Scybalocanthon sp. 01H")) #Select records of S. sexspilotus (sp 01H) from "morphospecies unification with SIB "DB-Col-2022-07-23"
    sc_sex <- mutate(sc_sex, scientificName1 = ifelse(scientificName1 == "Scybalocanthon sp. 01H", "Scybalocanthon sexspilotus", scientificName1)) # replace names
  records <- merge(records, sc_sex, by = intersect(names(records), names(sc_sex)), all = TRUE)# add records S. sexspilotus
  records <- records[!(records$scientificName1=="Sulcophanaeus auricollis" & records$coordinates == "-72.683333_5.433333"),] # remove atypical altitudinal records of S. auricollis
  records <- records[!(records$scientificName1=="Sulcophanaeus noctis" & records$coordinates %in% c("-75.55_4.716667", "-74.783333_10.95")),]
  records <- records[!(records$scientificName1=="Sulcophanaeus velutinus" & records$coordinates == "-78.10918_-1.212233"),] # remove atypical altitudinal records of S. velutinus
  records <- records[!(records$scientificName1 == "Sylvicanthon aequinoctialis" & records$coordinates == "-73.999611_10.895134"), ] # remove atypical records of S. aequinoctialis
  records <- records[!(records$scientificName1 == "Sylvicanthon aequinoctialis" & records$stateProvince %in% c("Amazonas", "Caquetá","Madre de Dios", "Nariño","Putumayo", "Vichada")), ] # remove S. aequinoctialis from Amazonas
  records <- records[!(records$scientificName1=="Sylvicanthon aequinoctialis" & records$decimalLongitude >= -72.9),] # remove records less than -72.9 Long S. aequinoctialis
  records <- records[!(records$scientificName1=="Sylvicanthon aequinoctialis" & records$decimalLatitude <= 0),] # remove records less than 0 Lat S. aequinoctialis
    s_proseni <- data.frame(subset(SIB_COL, scientificName1 == "Sylvicanthon proseni")) # Select records of Sylvicanthon proseni (sp 02H) from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, s_proseni, by = intersect(names(records), names(s_proseni)), all = TRUE) # add records of S proseni
  records <- records[!(records$scientificName1=="Sylvicanthon proseni" & records$coordinates == "-76.1_1.616666667"),]
    U_caucanus <- data.frame(subset(SIB_COL, scientificName1 == "Uroxys caucanus")) # Select records of U. caucanus from "morphospecies unification with SIB "DB-Col-2022-07-23"
  records <- merge(records, U_caucanus, by = intersect(names(records), names(U_caucanus)), all = TRUE) # add records of U caucanus
  records <- records[!(records$scientificName1 == "Uroxys caucanus" & records$stateProvince %in% c("Cundinamarca","Santander")), ] # remove U. caucanus from eastern cordillera
    u_pygma <- data.frame(subset(SIB_COL, scientificName1=="Uroxys pygmaeus")) # select records of U. pygmaeus from Entomological collections
  records <- merge(records, u_pygma, by=intersect(names(records), names(u_pygma)), all=TRUE) # add records of U. pygmaeus
  records <- records[!is.na(records$scientificName1),]
  
  saveRDS(records, file="registros_clean_2024-09.rds")
  write.table(records, "registros_clean_2024-09.txt")
#
# Morphospecies data: Entomological collections (IAvH-E, MEFLG, ICN, UPTC, CEUN-PSO)
  #
  msp_SIB <- data.frame(subset(SIB_COL, taxonRank == "Morfoespecie"))
  msp_project <- data.frame(subset(IBD_data, taxonRank == "Morphospecies IAvH"))
  msp_names <- sort(unique(msp_project$scientificName1))
  msp_records <- data.frame(msp_SIB[msp_SIB$scientificName1 %in% msp_names, ])
  msp_records <- msp_records[!(msp_records$scientificName1 =="Canthidium sp. 12H" & msp_records$coordinates =="-67.87083333_5.350833333"),]
  msp_records <- msp_records[!(msp_records$scientificName1=="Canthidium sp. 19H" & msp_records$coordinates <=-74),]
  msp_records <- msp_records[!(msp_records$scientificName1=="Uroxys gr. pauliani" & msp_records$coordinates == "-76.95691_3.69862"),]
  
  saveRDS(msp_records, file="msp_records.rds")
  write.table(msp_records, "msp_records.txt")
#
# merge dataset species/morphospecies from GBIF, project
#
  GBIF_final <- records[, c("scientificName1", "decimalLatitude", "decimalLongitude")]
  GBIF_final$source <- "GBIF"
  msp_final <- msp_records[, c("scientificName1", "decimalLatitude", "decimalLongitude")]
  msp_final$source <- "Collections"
  beta_db <- IBD_data[, c("scientificName1", "lat_all_points", "lon_all_points")]
  colnames(beta_db)[colnames(beta_db) == "lat_all_points"] <- "decimalLatitude"
  colnames(beta_db)[colnames(beta_db) == "lon_all_points"] <- "decimalLongitude"
  beta_db$source <- "Beta diversity"
  
  records_combined <- rbind(GBIF_final, msp_final, beta_db)
  
  saveRDS(records_combined, "./records_combined.rds")
  
  
################################################################################