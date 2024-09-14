setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata/GBIF_data")
# packages
  library(dplyr)
  library(ggplot2)
  library(rnaturalearth)
  library(rgdal)
  library(sf)
  library(readxl)
#
# GBIF dataset, morphospecies SIB-Col and upload records for species
    GBIF <- read_excel("GBIF_2024-02-13.xlsx", sheet = "GBIF_2024-02-13")
      GBIF$coordenadas <- paste(GBIF$decimalLongitude, GBIF$decimalLatitude, sep = "_")
    SIB_COL <- read_excel("DB-Col_2024-01-23.xlsx", sheet = "Conjunto-2020-07-23")
      SIB_COL$coordenadas <- paste(SIB_COL$decimalLongitude, SIB_COL$decimalLatitude, sep = "_")
      SIB_COL$scientificName1 <- paste(SIB_COL$identificationRemarks)
    at_aeneomicans <- read_excel("at_aeneomicans.xlsx", sheet = "at_aeneomicans")
      at_aeneomicans$coordenadas <- paste(at_aeneomicans$decimalLongitude, at_aeneomicans$decimalLatitude, sep = "_") 
      at_aeneomicans$scientificName1 <- paste(at_aeneomicans$genus, at_aeneomicans$specificEpithet)
      at_aeneomicans <- at_aeneomicans[!duplicated(at_aeneomicans$coordenadas), ]
    c_juvencus <- read_excel("c_juvencus.xlsx", sheet = "c_juvencus")
      c_juvencus$coordenadas <- paste(c_juvencus$decimalLongitude, c_juvencus$decimalLatitude, sep = "_")
      c_juvencus$scientificName1 <- paste(c_juvencus$scientificName)
      c_juvencus <- mutate(c_juvencus, scientificName1 = ifelse(scientificName1 == "Canthon juvencus Harold, 1868", "Canthon juvencus", scientificName1))
      c_juvencus <- mutate(c_juvencus, scientificName1 = ifelse(scientificName1 == "Canthon raripilus Bates, 1887", "Canthon juvencus", scientificName1))
      c_juvencus <- mutate(c_juvencus, scientificName1 = ifelse(scientificName1 == "BOLD:ABU7370", "Canthon juvencus", scientificName1))
      c_juvencus <- c_juvencus[!duplicated(c_juvencus$coordenadas), ]
    cps_corythus <- read_excel("c_corythus.xlsx", sheet = "c_corythus")
      cps_corythus$coordenadas <- paste(cps_corythus$decimalLongitude, cps_corythus$decimalLatitude, sep = "_")
      cps_corythus$scientificName1 <- paste(cps_corythus$genus, cps_corythus$specificEpithet)
      cps_corythus <- mutate(cps_corythus, scientificName1 = ifelse(scientificName1 == "Coprophanaeus telamon", "Coprophanaeus corythus", scientificName1)) # replace names
      cps_corythus <- cps_corythus[!duplicated(cps_corythus$coordenadas), ]
    d_guildingii <- read_excel("d_guildingii.xlsx", sheet = "d_guildingii")
      d_guildingii$coordenadas <- paste(d_guildingii$decimalLongitude, d_guildingii$decimalLatitude, sep = "_")
      d_guildingii$scientificName1 <- paste(d_guildingii$genus, d_guildingii$specificEpithet)
      d_guildingii <- d_guildingii[!duplicated(d_guildingii$coordenadas), ]
    di_agenor <- read_excel("d_agenor.xlsx", sheet = "d_agenor")
      di_agenor[, c("decimalLongitude", "decimalLatitude")] <- lapply(di_agenor[, c("decimalLongitude", "decimalLatitude")], as.numeric)
      di_agenor$coordenadas <- paste(di_agenor$decimalLongitude, di_agenor$decimalLatitude, sep = "_")
      di_agenor$scientificName1 <- paste(di_agenor$genus, di_agenor$specificEpithet)
      di_agenor <- di_agenor[!duplicated(di_agenor$coordenadas), ]
    d_gazella <- read_excel("d_gazella.xlsx", sheet = "d_gazella")
      d_gazella$coordenadas <- paste(d_gazella$decimalLongitude, d_gazella$decimalLatitude, sep = "_")
      d_gazella$scientificName1 <- paste(d_gazella$genus, d_gazella$specificEpithet)
      d_gazella <- mutate(d_gazella, scientificName1 = ifelse(scientificName1 == "Digitonthophagus NA", "Digitonthophagus gazella", scientificName1)) # replace NA for gazella in the specific epithet
      d_gazella <- d_gazella[!duplicated(d_gazella$coordenadas), ]
    eury_camero <- read_excel("eury_camero.xlsx", sheet = "Hoja1")
      eury_camero$coordenadas <- paste(eury_camero$decimalLongitude, eury_camero$decimalLatitude, sep = "_")
      eury_camero$scientificName1 <- paste(eury_camero$genus, eury_camero$specificEpithet)
    o_acuminatus <- read_excel("o_acuminatus.xlsx", sheet = "o_acuminatus")
      o_acuminatus$coordenadas <- paste(o_acuminatus$decimalLongitude, o_acuminatus$decimalLatitude, sep = "_")
      o_acuminatus$scientificName1 <- paste(o_acuminatus$genus, o_acuminatus$specificEpithet)
      o_acuminatus <- mutate(o_acuminatus, scientificName1 = ifelse(scientificName1 == "Onthophagus NA", "Onthophagus acuminatus", scientificName1)) # replace NA for acuminatus in the specific epithet
      o_acuminatus <- o_acuminatus[!duplicated(o_acuminatus$coordenadas), ]
      o_curvicornis <- read_excel("o_curvicornis.xlsx", sheet = "o_curvicornis")
    o_curvicornis$coordenadas <- paste(o_curvicornis$decimalLongitude, o_curvicornis$decimalLatitude, sep = "_")
      o_curvicornis$scientificName1 <- paste(o_curvicornis$genus, o_curvicornis$specificEpithet)
      o_curvicornis <- o_curvicornis[!duplicated(o_curvicornis$coordenadas), ]
    s_aequinoctialis <- read_excel("s_aequinoctialis.xlsx", sheet = "s_aequinoctialis")
      s_aequinoctialis$coordenadas <- paste(s_aequinoctialis$decimalLongitude, s_aequinoctialis$decimalLatitude, sep = "_")
      s_aequinoctialis$scientificName1 <- paste(s_aequinoctialis$genus, s_aequinoctialis$specificEpithet)
      s_aequinoctialis <- s_aequinoctialis[!duplicated(s_aequinoctialis$coordenadas), ]
#
# clean
    GBIF <- GBIF[GBIF$coordenadas != "-74.266667_4.05", ] # delete recrods of An. villosus in páramo
  registros <- merge(GBIF, at_aeneomicans, by = intersect(names(GBIF), names(at_aeneomicans)), all = TRUE) # add records of at_aeneomicans
  registros <- registros[registros$coordenadas != "-5.983333_-2.416667", ] # Remove records of A. irinus in the ocean
  registros <- registros[registros$coordenadas != "-74.855_2.7975", ] # Remove atypic records of A. irinus
    cdm_on <- data.frame(subset(SIB_COL, identificationRemarks == "Canthidium sp. 36H")) # Select records of Canthidium onitoides from "morphospecies unification with SIB "DB-Col-2022-07-23"
    cdm_on <- mutate(cdm_on, scientificName1 = ifelse(scientificName1 == "Canthidium sp. 36H", "Canthidium onitoides", scientificName1))  # replace names
    cdm_on <- cdm_on[!duplicated(cdm_on$coordenadas), ]
  registros <- merge(registros, cdm_on, by = intersect(names(registros), names(cdm_on)), all = TRUE) # add records for C onitoides
    ctn_fu <- subset(SIB_COL, identificationRemarks =="Canthon fulgidus")# Select records of Canthon fulgidus from "morphospecies unification with SIB "DB-Col-2022-07-23"
    ctn_fu <- mutate(ctn_fu, scientificName1 = ifelse(scientificName1 == "Canthon fulgidus", "Canthon fulgidus martinezi", scientificName1)) # replace names
    ctn_fu <- ctn_fu[!duplicated(ctn_fu$coordenadas), ]
  registros <- merge(registros, ctn_fu, by = intersect(names(registros), names(ctn_fu)), all=TRUE)# add records for C. fulgidus martinezi
  registros <- merge(registros, c_juvencus, by = intersect(names(registros), names(c_juvencus)), all = TRUE) # add records for C. juvencus
  registros <- registros[!(registros$scientificName1=="Canthon juvencus" & registros$coordenadas %in% c("-74.25_4.05", "-75.078596_6.913269", "-72.536611_11.145056")), ] # Remove atypical altitudinal records of C. juvencus
    ctn_li <- data.frame(subset(SIB_COL, scientificName1 == "Canthon lituratus")) # Select records of Canthon lituratus from "morphospecies unification with SIB "DB-Col-2022-07-23"
    ctn_li <- ctn_li[!duplicated(ctn_li$coordenadas), ]
  registros <- merge(registros, ctn_li, by = intersect(names(registros), names(ctn_li)), all = TRUE) # add records for C lituratus
  registros <- registros[!(registros$scientificName1=="Canthon lituratus" & registros$coordenadas %in% c("-75078_10823","-74.84_9741", "-59.83333333_-2.416666667", "-74.083333_11.116667")),] # Remove atypical altitudinal records of C. lituratus
  registros <- registros[!(registros$scientificName1=="Canthon luteicollis" & registros$stateProvince %in% c("Arauca", "Casanare", "Vichada")),] # Remove records of C. luteicollis from Orinoquia
  registros <- registros[!(registros$scientificName1=="Canthon luteicollis" & registros$coordenadas == "-77.6_-0.46667"),] # Remove atypical altitudinal records of de C. luteicollis
  registros <- registros[!(registros$scientificName1=="Canthon pallidus" & registros$coordenadas =="-74.166667_2.666667"),]
  registros <- registros[!(registros$scientificName1=="Canthon semiopacus" & registros$coordenadas =="-79.24529_-3.99067"),]
  registros <- registros[!(registros$scientificName1=="Canthon septemmaculatus" & registros$decimalLatitude <= 5),] # Remove records than 5 Longitude for C. septemmaculatus
    ctn_su <- data.frame(subset(SIB_COL, identificationRemarks == "Canthon subhyalinus")) # Select records of Canthon subhyalinus from "morphospecies unification with SIB "DB-Col-2022-07-23"
    ctn_su <- ctn_su[!duplicated(ctn_su$coordenadas), ]
  registros <- merge(registros, ctn_su, by = intersect(names(registros), names(ctn_su)), all = TRUE) # add records for Canthon subhyalinus
  registros <- registros[registros$coordenadas != "-74831_4936", ] # Remove atypical records of C. subhyalinus
  registros <- registros[!(registros$scientificName1=="Canthon subhyalinus" & registros$coordenadas =="-75.07859583_6.913268902"),] # Remove atypical altitudinal records of C. subhyalinus
  registros <- registros[!(registros$scientificName1=="Canthon subhyalinus" & registros$decimalLatitude <= 1.7),] # Remove records less than 1.7 Longitude of C. subhyalinus
  registros <- merge(registros, cps_corythus, by = intersect(names(registros), names(cps_corythus)), all = TRUE) # add records of Coporphanaeus corythus (with synonims and changes in taxonomy)
  registros <- registros[!(registros$scientificName1=="Coprophanaeus corythus" & registros$stateProvince %in% c("Meta", "Casanare", "Vichada")),] # remove records of C. corythus from Orinoquia
  registros <- registros[!(registros$scientificName1=="Coprophanaeus corythus" & registros$decimalLatitude >= 23),] # delete records less than 23 Latitude of C. corythus
  registros <- registros[!(registros$scientificName1 == "Coprophanaeus jasius" & registros$stateProvince %in% c("Arauca", "Antioquia", "Canindeyu", "Sucre", "Atlántico")), ] # remove records of C. jasius from Atlantico department
  registros <- registros %>% mutate(scientificName1 = ifelse(scientificName1 == "Coprophanaeus ohausi" & !is.na(stateProvince) & stateProvince == "Antioquia", "Coprophanaeus morenoi", scientificName1)) # Change ID speices from C. ohausi (Mistake ID) to C. morenoi
  registros <- registros[!(registros$scientificName1=="Coprophanaeus morenoi" & registros$coordenadas == "-78.42578_-1.39315"),] # Remove atypical altitudinal records of C. morenoi
  registros <- registros[!(registros$scientificName1 == "Coprophanaeus telamon" &! is.na(registros$stateProvince) & registros$stateProvince %in% c("Antioquia","Caldas", "Cesar", "Córdoba", "Cundinamarca","La Guajira", "Santander", "Sucre", "Valle del Cauca")), ] # remove some departments from the distribution of C. telamon
  registros <- registros[!(registros$scientificName1 == "Coprophanaeus telamon" &! is.na(registros$county) & registros$county == "Puerto Boyacá"), ] # Remove C. telamon from Magdalena Valley (It really is C. corythys)
  registros <- registros[!(registros$scientificName1 == "Deltochilum carinatum" &! is.na(registros$country) & registros$country == "Costa Rica"), ] # Remove D. carinatum from Costa Rica
  registros <- registros[!(registros$scientificName1 == "Deltochilum carinatum" & registros$coordenadas %in% c("-76.1_1.333333", "-76.544889_1.243611")), ] # Remove D. carinatum from Andes
  registros <- merge(registros, d_guildingii, by = intersect(names(registros), names(d_guildingii)), all = TRUE) # add records for D. guildingii
  registros <- registros[!(registros$scientificName1=="Deltochilum guildingii" & registros$stateProvince <= "CaquetÃ¡"),] # remove D. guildingii from Caquetá
  registros <- registros[!(registros$scientificName1=="Deltochilum guildingii" & registros$coordenadas=="-74.083333_11.116667"), ] # remove atypical altitudinal records of D. guildingii
  registros <- registros[!(registros$scientificName1=="Deltochilum hypponum" & registros$coordenadas %in% c("-75.428647_0.838333", "-75.246667_2.935278", "-76.633333_1.133333", "-73.388056_4.5925", "-73.408611_5.75")), ] # remove atypical altitudinal records of D. hypponum
    del_mex <- data.frame(subset(SIB_COL, identificationRemarks == "Deltochilum mexicanum")) # Select records of Deltochilum mexicanum from "morphospecies unification with SIB "DB-Col-2022-07-23"
    del_mex <- del_mex[!duplicated(del_mex$coordenadas), ]
  registros <- merge(registros, del_mex, by = intersect(names(registros), names(del_mex)), all = TRUE) # add records for D. mexicanum
  registros <- registros[!(registros$scientificName1=="Deltochilum mexicanum" & registros$coordenadas %in% c("-71.8375_6.513333333", "-73.38805556_4.5925", "-78.25_1.283333333", "-74.72043_2.71329")), ] # remove atypical altitudinal records of D. burmeisteri
  registros <- registros[!(registros$scientificName1=="Deltochilum mexicanum" & !is.na(registros$country) & registros$country %in% c("Guatemala","Mexico")), ] # Remove D. burmeisteri from Guatemala and Mexico 
  registros <- registros[!(registros$scientificName1=="Deltochilum molanoi" & registros$coordenadas == "-76.633333_2.283333"),] # remove points greater than 76 Longitude D, molanoi
  registros <- registros[!(registros$scientificName1=="Deltochilum orbiculare" & registros$decimalLongitude <= -76.4),] # remove points greater than 76 Longitude D. orbiculare
  registros <- merge(registros, di_agenor, by = intersect(names(registros), names(di_agenor)), all = TRUE) # add records of D. agenor
  registros <- registros[!(registros$scientificName1=="Dichotomius agenor" & registros$decimalLongitude <= -81),] # remove points greater than  81 Long D agenor
  registros <- registros[!(registros$scientificName1=="Dichotomius agenor" &! is.na(registros$stateProvince) & registros$stateProvince %in% c("Meta","Casanare", "Vichada", "BoyacÃ¡")), ] # Remove D. agenor departments based on Montoya et al. 2021
  registros <- registros[!(registros$scientificName1=="Dichotomius agenor" & registros$coordenadas %in% c("-74.84_9741", "-74907_9027", "-74831_4936", "NA_NA", "-74907_9027","-75078_10823", "-74831_4936", "-74.06785_11.0973")), ] # remove atypical records of D. agenor
    di_bel <- data.frame(subset(SIB_COL, identificationRemarks == "Dichotomius belus")) # Select records of Dichotomius belus from "morphospecies unification with SIB "DB-Col-2022-07-23"
    di_bel <- di_bel[!duplicated(di_bel$coordenadas), ]
  registros <- merge(registros, di_bel, by = intersect(names(registros), names(di_bel)), all = TRUE) # add records for D. belus
  registros <- registros[!(registros$scientificName1=="Dichotomius belus" & registros$decimalLatitude == 4936),] # remove atypical records of D. belus
  registros <- registros[!(registros$scientificName1=="Dichotomius belus" & registros$coordenadas %in% c("-74.53146_6.27808", "-74.73175_5.73692", "-74.554144_6.460248", "-74.572882_6.494321")), ] # remove atypical records of  D. belus
  registros <- registros[!(registros$scientificName1=="Dichotomius belus" & registros$decimalLongitude >= -72),] # remove points greater than -72 Long D. belus
  registros <- registros[!(registros$scientificName1 == "Dichotomius belus" &! is.na(registros$stateProvince) & registros$stateProvince %in% c("Cesar", "Atlántico", "Magdalena","La Guajira")), ] # Remove D. belus departments
    di_gan <- data.frame(subset(SIB_COL, identificationRemarks == "Dichotomius sp. 10H")) # Select records of Dichotomius gandinii from "morphospecies unification with SIB "DB-Col-2022-07-23"
    di_gan <- mutate(di_gan, scientificName1 = ifelse(scientificName1 == "Dichotomius sp. 10H", "Dichotomius gandinii", scientificName1)) # replace names
    di_gan <- di_gan[!duplicated(di_gan$coordenadas), ]
  registros <- merge(registros, di_gan, by = intersect(names(registros), names(di_gan)), all = TRUE) # add records of D. gandinii
    di_glo <- data.frame(subset(SIB_COL, identificationRemarks == "Dichotomius globulus")) # Select records of Dichotomius globulus from "morphospecies unification with SIB "DB-Col-2022-07-23"
    di_glo <- di_glo[!duplicated(di_glo$coordenadas), ]
  registros <- merge(registros, di_glo, by = intersect(names(registros), names(di_glo)), all = TRUE) # add records of D. globulus
    di_mam <- data.frame(subset(SIB_COL, identificationRemarks == "Dichotomius mamillatus")) # Select records of Dichotomius mamillatus from "morphospecies unification with SIB "DB-Col-2022-07-23"
    di_mam <- di_mam[!duplicated(di_mam$coordenadas), ]
  registros <- merge(registros, di_mam, by = intersect(names(registros), names(di_mam)), all = TRUE) # add records of D mamillatus
  registros <- registros[!(registros$scientificName1=="Dichotomius nisus" & registros$country == "Brazil"),] # Remove D nisus from Basil, because we don't like include Amazonia. D nisus inhabit open areas such as Orinoquia, el Cerrado, Catinga ...
    di_pro <- data.frame(subset(SIB_COL, identificationRemarks == "Dichotomius protectus")) #Select records of Dichotomius protcetus from "morphospecies unification with SIB "DB-Col-2022-07-23"
    di_pro <- di_pro[!duplicated(di_pro$coordenadas), ]
  registros <- merge(registros, di_pro, by = intersect(names(registros), names(di_pro)), all = TRUE) # add records of D protectus
  registros <- registros[!(registros$scientificName1 =="Dichotomius quadrilobatus" & registros$coordenadas =="-77_1.15"),]
  registros <- registros[!(registros$scientificName1 == "Dichotomius quinquelobatus" &! is.na(registros$stateProvince) & registros$stateProvince %in% c("Cauca", "Valle del Cauca", "Quindío","Risaralda", "Caldas")), ] # remove D. quinquelobatus departments
  registros <- registros[!(registros$scientificName1=="Dichotomius quinquelobatus" & registros$coordenadas %in% c("-78.433692_-0.367544", "-79.24529_-3.99067", "-77.21667_0.05")),] # remove atypical points of D. quinquelobatus
    di_val <- data.frame(subset(SIB_COL, identificationRemarks == "Dichotomius validipilosus")) # Select records of Dichotomius validipilosus from "morphospecies unification with SIB "DB-Col-2022-07-23"
  registros <- merge(registros, di_val, by = intersect(names(registros), names(di_val)), all = TRUE) # add records of D validipilosus
  registros <- registros[!(registros$scientificName1=="Dichotomius validipilosus" & registros$decimalLongitude >= -60),] # Remove records less than -60 Long D. validipilosus
  registros <- merge(registros, d_gazella, by = intersect(names(registros), names(d_gazella)), all = TRUE) # add records of D gazella GBIF
  registros <- registros[!(registros$scientificName1=="Digitonthophagus gazella" & registros$country != "Colombia"), ] # records from Colombia
  registros <- registros[!(registros$scientificName1=="Digitonthophagus gazella" & registros$coordenadas %in% c("-72.536611_11.145056", "-76.068569_3.105836")),] # remove atypical records of D. gazella
    eu_car <- data.frame(subset(eury_camero, scientificName == "Eurysternus caribaeus")) # select records of E. caribaeus Camero & Lobo 2012
    eu_car <- eu_car[!duplicated(eu_car$coordenadas), ]
  registros <- merge(registros, eu_car, by = intersect(names(registros), names(eu_car)), all = TRUE) # add records E. caribaeus Camero & Lobo 2012
  registros <- registros[!(registros$scientificName1=="Eurysternus caribaeus" & registros$coordenadas %in% c("-52.95_5.36666666666667", "-75.7166666666667_4.11666666666667", "-72.0833333333333_7.95566666666667")),] # remove atypical altitudinal records of E. caribaeus
    eu_cay <- data.frame(subset(eury_camero, scientificName == "Eurysternus cayennensis")) # select records of E. cayennensis Camero & Lobo 2012
    eu_cay <- eu_cay[!duplicated(eu_cay$coordenadas), ]
  registros <- merge(registros, eu_cay, by = intersect(names(registros), names(eu_cay)), all = TRUE) #  add records of E. cayennesis
  registros <- registros[!(registros$scientificName1=="Eurysternus cayennensis" & registros$coordenadas %in% c("-78.5666666666667_-0.0588333333333333", "-60.9166666666667_5.0505", "-52.95_5.36666666666667")), ] # remove atypical altitudinal records of E. cayennensis
  registros <- registros[!(registros$scientificName1=="Eurysternus cayennensis" &! is.na(registros$stateProvince) & registros$stateProvince %in% c("BOLÍVAR","CHOCÓ")), ] #  Remove E. cayennesis from departments of Caribe
  registros <- registros[!(registros$scientificName1=="Eurysternus cayennensis" & registros$decimalLatitude <= -20),] # Remove records less tan -20 Lat E. cayennensis
  registros <- registros[!(registros$scientificName1=="Eurysternus contractus" & registros$coordenadas %in% c("-77.3333_6.00667", "-70.533333_0.233333", "-78.936729_-1.07175", "-79.005517_-2.900755", "-75.666553_1.508061", "-77.216667_0.05")), ] # remove atypical records of E. contractus
    eu_foe <- data.frame(subset(eury_camero, scientificName == "Eurysternus foedus")) # select records of E. foedus Camero & Lobo 2012
    eu_foe <- eu_foe[!duplicated(eu_foe$coordenadas),]
  registros <- merge(registros, eu_foe, by = intersect(names(registros), names(eu_foe)), all = TRUE) # add records of E. foedus
  registros <- registros[!(registros$scientificName1=="Eurysternus foedus" & ! is.na(registros$locality) & registros$locality %in% c("Vacas", "Chugchilán", "Río Frío, Parque Nacional El Tamá")), ] # remove atypical altitudinal records of E. foedus
  registros <- registros[!(registros$scientificName1=="Eurysternus foedus" & registros$decimalLatitude <= -20),] # remove records less than -20 Lat E. foedus
  registros <- registros[!(registros$scientificName1=="Eurysternus foedus" & registros$coordenadas == "-67.3333333333333_-16.4666666666667"),] # remove atypical altitudinal records of E. foedus
    eu_hyp <- data.frame(subset(eury_camero, scientificName == "Eurysternus hypocrita")) # select records of E. hypocrita Camero & Lobo 2012
    eu_hyp <- eu_hyp[!duplicated(eu_hyp$coordenadas),]
  registros <- merge(registros, eu_hyp, by = intersect(names(registros), names(eu_hyp)), all = TRUE) # add records of E. hypocrita
  registros <- registros[!(registros$scientificName1=="Eurysternus hypocrita" & registros$coordenadas %in% c("-75.5666666666667_-10.55", "-61.4333333333333_5.98333333333333", "-52.95_5.36666666666667")), ] # remove atypical altitudinal records of E. hypocrita
  registros <- registros[!(registros$scientificName1 == "Eurysternus hypocrita" &! is.na(registros$stateProvince) & registros$stateProvince %in% c("Antioquia","Santander")), ] # Remove E. hypocrita from Magdalena Valley
  registros <- registros[!(registros$scientificName1=="Eurysternus hypocrita" & registros$decimalLatitude <= -20),] #remove records less than -20 Lat E. hypocrita
    eu_mar <- data.frame(subset(eury_camero, scientificName == "Eurysternus marmoreus")) # select records of E. marmoreus Camero & Lobo 2012
    eu_mar <- eu_mar[!duplicated(eu_mar$coordenadas),]
  registros <- merge(registros, eu_mar, by = intersect(names(registros), names(eu_mar)), all = TRUE) # add records of E. marmoreus
  registros <- registros[!(registros$scientificName1=="Eurysternus marmoreus" & registros$coordenadas %in% c("-75.4333333333333_4.75166666666667", "-73.0833333333333_5.91666666666667","-70.9833333333333_-13.4833333333333","-78.1166666666667_1.66666666666667", "-74.033333_11.333333")), ] # remove atypical altitudinal points of E. marmoreus
    eu_mex <- data.frame(subset(eury_camero, scientificName1 == "Eurysternus mexicanus")) # select records of E. mexicanus Camero & Lobo 2012
    eu_mex <- eu_mex[!duplicated(eu_mex$coordenadas),]
  registros <- merge(registros, eu_mex, by = intersect(names(registros), names(eu_mex)), all = TRUE) # add records of E. mexicanus
  registros <- registros[!(registros$scientificName1=="Eurysternus mexicanus" & registros$coordenadas %in% c("-74.1166666666667_4.46666666666667", "-71.9833333333333_7.95166666666667")), ]# remove atypical altitudinal records of E. mexicanus
    eu_ple <- data.frame(subset(eury_camero, scientificName == "Eurysternus plebejus")) # select records of E. plebejus Camero & Lobo 2012
    eu_ple <- eu_ple[!duplicated(eu_ple$coordenadas),]
  registros <- merge(registros, eu_ple, by = intersect(names(registros), names(eu_ple)), all = TRUE) # add records of E. plebejus
  registros <- registros[!(registros$scientificName1=="Eurysternus plebejus" & registros$coordenadas %in% c("-75.5166666666667_6.78333333333333", "-79.195752_-4.004985", "-79.199326_-3.986796", "-65.5833333333333_-17.5166666666667")), ] # remove atypical altitudinal records of E. plebejus
    eu_sqa <- data.frame(subset(eury_camero, scientificName == "Eurysternus squamosus")) # select records of E. squamosus Camero & Lobo 2012
  registros <- merge(registros, eu_sqa, by = intersect(names(registros), names(eu_sqa)), all = TRUE) # add records of E. squamosus
  registros <- registros[!(registros$scientificName1=="Eurysternus squamosus" & registros$coordenadas %in% c("-77.6166666666667_0.816666666666667", "-76.1_1.333333")),] #elimina puntos altitudinales atipicos de E. squamosus
    eu_wit <- data.frame(subset(eury_camero, scientificName == "Eurysternus wittmerorum")) # select records of E. wittmerorum Camero & Lobo 2012
    eu_wit <- eu_wit[!duplicated(eu_wit$coordenadas),]
  registros <- merge(registros, eu_wit, by = intersect(names(registros), names(eu_wit)), all = TRUE) # add records of E. wittmerorum
  registros <- registros[!(registros$scientificName1=="Eurysternus wittmerorum" & registros$coordenadas %in% c("-75.55_4.716667","-79.199326_-3.986796")),]
  registros <- registros[!(registros$scientificName1=="Eurysternus wittmerorum" & registros$decimalLatitude >= 4.8),] # Remove E. wittmerorum from Boyaca
  registros <- registros[!(registros$scientificName1 == "Gromphas aeruginosa" & registros$coordenadas %in% c("-71.46_7.57", "-6.333333_-12.833333")), ] # remove atypical records of G. aeruginosa follow Cupello et al. 2013
  registros <- registros[!(registros$scientificName1=="Gromphas lemoinei" & registros$coordenadas == "-74.081944_4.61"),] # Remove records of G aeruginosa from catinga ...
    on_bre <- data.frame(subset(SIB_COL, identificationRemarks == "Ontherus brevicollis")) # Select records of Ontherus brevicollis from "morphospecies unification with SIB "DB-Col-2022-07-23"
    on_bre <- on_bre[!duplicated(on_bre$coordenadas), ]
  registros <- merge(registros, on_bre, by = intersect(names(registros), names(on_bre)), all = TRUE) # add records of O. brevicollis
  registros <- registros[!(registros$scientificName1=="Ontherus brevicollis" & registros$coordenadas %in% c("-75.246667_2.935278", "-73.38805556_4.5925", "-73.3808_6.93708", "-73.38349_4.56448")), ] #Eliminar puntos de altura atipicos O. brevicollis
  registros <- registros[!(registros$scientificName1=="Ontherus diabolicus" & registros$coordenadas == "-72.683333_5.433333"),] # remove atypical altitudinal records of O. diabolicus
    on_inc <- data.frame(subset(SIB_COL, identificationRemarks == "Ontherus incisus")) # Select records of Ontherus incisus from "morphospecies unification with SIB "DB-Col-2022-07-23"
    on_inc <- on_inc[!duplicated(on_inc$coordenadas), ]
  registros <- merge(registros, on_inc, by = intersect(names(registros), names(on_inc)), all = TRUE) # add records of O. incisus
  registros <- registros[!(registros$scientificName1=="Ontherus incisus" & registros$stateProvince == "Santander"),] # remove O. incicus from Santander
  registros <- registros[!(registros$scientificName1=="Ontherus kirschii" & registros$coordenadas %in% c("-72.840833_7.438889", 	"-73.388056_4.5925")),] # remove atypical altitudinal records of O. kirschii
    on_lun <- data.frame(subset(SIB_COL, identificationRemarks == "Ontherus lunicollis")) # Select records of Ontherus lunicollis from "morphospecies unification with SIB "DB-Col-2022-07-23"
    on_lun <- on_lun[!duplicated(on_lun$coordenadas), ]
  registros <- merge(registros, on_lun, by = intersect(names(registros), names(on_lun)), all = TRUE) #a add records of O. lunicollis
  registros <- registros[!(registros$scientificName1=="Ontherus lunicollis" & registros$coordenadas %in% c("-75.246667_2.935278", "-72.61869444_6.66511111", "-75.495361_4.469028")), ] # remove atypical altitudinal records of O. lunicollis
  registros <- registros[!(registros$scientificName1=="Ontherus pubens" & registros$coordenadas == "-78.42578_-1.39315"),] # remove atypical altitudinal records of O. pubens
  registros <- merge(registros, o_acuminatus, by = intersect(names(registros), names(o_acuminatus)), all = TRUE) # add records of O. acuminatus
  registros <- registros[!(registros$scientificName1=="Onthophagus acuminatus" & registros$decimalLongitude >= -72.7),] # remove records less than -73 Lat O acuminatus
  registros <- registros[!(registros$scientificName1=="Onthophagus acuminatus" & registros$coordenadas %in% c("-73.1_9.85", "-73.458333_5.746667", "-76.068569_3.105836")),] #Eliminar puntos de altura atipicos en O. acuminatus
    o_bidentatus <- data.frame(subset(SIB_COL, identificationRemarks == "Onthophagus bidentatus")) # Select records of Onthophagus bidentatus from "morphospecies unification with SIB "DB-Col-2022-07-23"
    o_bidentatus <- o_bidentatus[!duplicated(o_bidentatus$coordenadas), ]
  registros <- merge(registros, o_bidentatus, by = intersect(names(registros), names(o_bidentatus)), all = TRUE) # add records of O. bidentatus
  registros <- registros[!(registros$scientificName1 == "Onthophagus bidentatus" & registros$coordenadas %in% c("-74907_9027", "-74.84_9741", "-75078_10823")), ] # remove atypical altitudinal records of O. bidendatus
  registros <- merge(registros, o_curvicornis, by = intersect(names(registros), names(o_curvicornis)), all = TRUE) # add records of O. curvicornis
  registros <- registros[!(registros$scientificName1=="Onthophagus curvicornis" & registros$decimalLatitude >= 10),] # remove records greater than 10 Lat O curvicornis
  registros <- registros[!(registros$scientificName1=="Onthophagus curvicornis" & registros$decimalLongitude >= -70),] # remove records gerater than 70 Long O curvicornis
  registros <- registros[!(registros$scientificName1=="Onthophagus curvicornis" & registros$coordenadas %in% c("-72_4", "-71.86113_6.47894", "-72.309928_5.608447", "-72.987472_5.946583")),] # remove atypical altitudinal records of O. curvicornis
    o_marg <- data.frame(subset(SIB_COL, identificationRemarks == "Onthophagus marginicollis")) #s Select records of Onthophagus marginicollis from "morphospecies unification with SIB "DB-Col-2022-07-23"
    o_marg <- o_marg[!duplicated(o_marg$coordenadas), ]
  registros <- merge(registros, o_marg, by = intersect(names(registros), names(o_marg)), all = TRUE) # add records of O. marginicollis
  registros <- registros[!(registros$scientificName1 == "Onthophagus marginicollis" & registros$decimalLatitude %in% c("9027", "9741")), ] # remove atypical records of O. marginicollis
    o_tran <- data.frame(subset(SIB_COL, identificationRemarks == "Onthophagus transisthmius")) # Select records of Onthophagus transisthmius from "morphospecies unification with SIB "DB-Col-2022-07-23"
    o_tran <- o_tran[!duplicated(o_tran$coordenadas), ]
  registros <- merge(registros, o_tran, by = intersect(names(registros), names(o_tran)), all = TRUE) # Add records of O. transithmius
    ox_ebe <- data.frame(subset(SIB_COL, identificationRemarks == "Oxysternon ebeninum")) # Select records of Oxysternon ebeninum from "morphospecies unification with SIB "DB-Col-2022-07-23"
    ox_ebe <- ox_ebe[!duplicated(ox_ebe$coordenadas), ]
  registros <- merge(registros, ox_ebe, by = intersect(names(registros), names(ox_ebe)), all = TRUE) # add records of O. ebeninum
    ox_sil <- data.frame(subset(SIB_COL, identificationRemarks == "Oxysternon silenus")) # Select records of Oxysternon silenus from "morphospecies unification with SIB "DB-Col-2022-07-23"
    ox_sil <- ox_sil[!duplicated(ox_sil$coordenadas), ]
  registros <- merge(registros, ox_sil, by = intersect(names(registros), names(ox_sil)), all = TRUE) # add records of O. silenus
  registros <- registros[!(registros$scientificName1=="Phanaeus bispinus" & registros$coordenadas=="-79.199939_-3.982992"),]
  registros <- registros[!(registros$scientificName1 == "Phanaeus chalcomelas" & registros$coordenadas %in% c("-78.433692_-0.367544", "-79.24529_-3.99067")), ] # remove atypical altitudinal records P. chalcomelas
  registros <- registros[!(registros$scientificName1=="Phanaeus haroldi" & registros$decimalLongitude >= -65),] # remove records less than -60 Long of P. haroldi
  registros <- registros[!(registros$scientificName1=="Phanaeus haroldi" & registros$coordenadas %in% c("-72.691667_5.434722", "-72.691667_5.434722")),] # remove atypical altitudinal records of P. haroldi
  registros <- registros[!(registros$scientificName1=="Phanaeus meleagris" & registros$coordenadas =="-79.199326_-3.986796"),]
  registros <- registros[!(registros$scientificName1 == "Phanaeus pyrois" & registros$country != "Colombia"), ] #extraer solo información de P. pyrois para Colombia
  registros <- registros[!(registros$scientificName1 == "Pseudocanthon perplexus"),] # remove records from "registros" of P. perplexus, and we work with unification with SIB "DB-Col-2022-07-23" 
    ps_per <- data.frame(subset(SIB_COL, identificationRemarks == "Pseudocanthon sp. 01H")) # Select records of Pseudocanthon perplexus (sp 01H) from "morphospecies unification with SIB "DB-Col-2022-07-23"
    ps_per <- mutate(ps_per, scientificName1 = ifelse(scientificName1 == "Pseudocanthon sp. 01H", "Pseudocanthon perplexus", scientificName1)) # replace names
    ps_per <- ps_per[!duplicated(ps_per$coordenadas), ]
  registros <- merge(registros, ps_per, by = intersect(names(registros), names(ps_per)), all = TRUE)# add records of P. perplexus
  registros <- registros[!(registros$scientificName1== "Pseudocanthon perplexus" & registros$decimalLatitude== "9027"),] # remove atypical altitudinal records of Pseudocanthon perplexus
    ps_xan <- data.frame(subset(SIB_COL, identificationRemarks == "Pseudocanthon sp. 02H")) # Select records of Pseudocanthon xanthurus (sp 02H) from "morphospecies unification with SIB "DB-Col-2022-07-23"
    ps_xan <- mutate(ps_xan, scientificName1 = ifelse(scientificName1 == "Pseudocanthon sp. 02H", "Pseudocanthon xanthurus", scientificName1)) # replace names
    ps_xan <- ps_xan[!duplicated(ps_xan$coordenadas), ]
  registros <- merge(registros, ps_xan, by = intersect(names(registros), names(ps_xan)), all = TRUE)# add records of P. xanthurus
  registros <- registros[!(registros$scientificName1=="Scatimus ovatus" & registros$decimalLongitude >= -72.4),] # remove records less than -72.4 Long S. ovatus
  registros <- registros[!(registros$scientificName1=="Scatimus ovatus" & registros$decimalLatitude <= 4.9),] # remove records less than 4.9 Lat S. ovatus
  registros <- registros[!(registros$scientificName1=="Scybalocanthon kelleri" & registros$coordenadas=="-73.595_3.045556"), ]# remove atypical altitudinal records of S. kelleri
  registros <- registros[!(registros$scientificName1=="Scybalocanthon darlingtoni" & registros$decimalLatitude >= 11),] # remove records greater than 11 Lat S. darlingtoni
    sc_dar <- data.frame(subset(SIB_COL, scientificName == "Scybalocanthon darlingtoni"))# Select records of Scybalocanton darlingtoni from "morphospecies unification with SIB "DB-Col-2022-07-23"
  registros <- merge(registros, sc_dar, by = intersect(names(registros), names(sc_dar)), all = TRUE)# add records of S. sexspilotus
    sc_mar <- data.frame(subset(SIB_COL, scientificName == "Scybalocanthon martinezi"))# Select records of S. martinezi from "morphospecies unification with SIB "DB-Col-2022-07-23"
    sc_mar <- sc_mar[!duplicated(sc_mar$coordenadas), ]
  registros <- merge(registros, sc_mar, by = intersect(names(registros), names(sc_mar)), all = TRUE)# add recrods of S. martinezi
    sc_sex <- data.frame(subset(SIB_COL, identificationRemarks == "Scybalocanthon sp. 01H")) #Select records of S. sexspilotus (sp 01H) from "morphospecies unification with SIB "DB-Col-2022-07-23"
    sc_sex <- mutate(sc_sex, scientificName1 = ifelse(scientificName1 == "Scybalocanthon sp. 01H", "Scybalocanthon sexspilotus", scientificName1)) # replace names
    sc_sex <- sc_sex[!duplicated(sc_sex$coordenadas), ]
  registros <- merge(registros, sc_sex, by = intersect(names(registros), names(sc_sex)), all = TRUE)# add records S. sexspilotus
  registros <- registros[!(registros$scientificName1=="Sulcophanaeus auricolis" & registros$coordenadas == "-72.683333_5.43333"),] # remove atypical altitudinal records of S. auricollis
  registros <- registros[!(registros$scientificName1=="Sulcophanaeus noctis" & registros$coordenadas =="-75.55_4.716667"),]
  registros <- registros[!(registros$scientificName1=="Sulcophanaeus velutinus" & registros$coordenadas == "-78.10918_-1.212233"),] # remove atypical altitudinal records of S. velutinus
  registros <- registros[registros$scientificName1 != "Sylvicanthon aequinoctialis", ]
  registros <- merge(registros, s_aequinoctialis, by = intersect(names(registros), names(s_aequinoctialis)), all = TRUE) # add records ofS aequinoctialis
  registros <- registros[!(registros$scientificName1 == "Sylvicanthon aequinoctialis" & registros$coordenadas %in% c("-74831_4936", "-74.84_9741", "-73.999611_10.895134")), ] # remove atypical records of S. aequinoctialis
  registros <- registros[!(registros$scientificName1 == "Sylvicanthon aequinoctialis" &! is.na(registros$stateProvince) & registros$stateProvince %in% c("Amazonas", "CaquetÃ¡","Madre de Dios", "NariÃ±o","Putumayo", "Vichada")), ] # remove S. aequinoctialis from Amazonas
  registros <- registros[!(registros$scientificName1=="Sylvicanthon aequinoctialis" & registros$decimalLongitude >= -72.9),] # remove records less than -72.9 Long S. aequinoctialis
    s_proseni <- data.frame(subset(SIB_COL, scientificName1 == "Sylvicanthon proseni")) # Select records of Sylvicanthon proseni (sp 02H) from "morphospecies unification with SIB "DB-Col-2022-07-23"
    s_proseni <- s_proseni[!duplicated(s_proseni$coordenadas), ]
  registros <- merge(registros, s_proseni, by = intersect(names(registros), names(s_proseni)), all = TRUE) # add records of S proseni
  registros <- registros[!(registros$scientificName1=="Sylvicanthon proseni" & registros$coordenadas == "-76.1_1.616666667"),]
    U_caucanus <- data.frame(subset(SIB_COL, scientificName1 == "Uroxys caucanus")) # Select records of U. caucanus from "morphospecies unification with SIB "DB-Col-2022-07-23"
    U_caucanus <- U_caucanus[!duplicated(U_caucanus$coordenadas), ]
  registros <- merge(registros, U_caucanus, by = intersect(names(registros), names(U_caucanus)), all = TRUE) # add records of U caucanus
  registros <- registros[!(registros$scientificName1 == "Uroxys caucanus" &! is.na(registros$stateProvince) & registros$stateProvince %in% c("Cundinamarca","Santander")), ] # remove U. caucanus from eastern cordillera
  registros <- registros[!is.na(registros$scientificName1),]
  
#saveRDS(registros, file="registros.rds")
#
# Check records by plot #
  south_america <- ne_countries(scale = "medium", continent = "South America", returnclass = "sf")
  central_america <- ne_countries(scale = "medium", continent = "North America", returnclass = "sf")
  central_america <- central_america[central_america$name %in% c("Belize", "Costa Rica", "El Salvador", "Mexico", "Guatemala", "Honduras", "Nicaragua", "Panama"), ]
  col <- ne_states(country = "colombia", returnclass = "sf")
  # subset species from GBIF
    sp <- data.frame(subset(registros,  scientificName1 == "Canthon fulgidus martinezi"))
    coordinates <- sp[,c("decimalLongitude", "decimalLatitude")]
    coordinates <- na.omit(coordinates)
    coordinates <- coordinates[!duplicated(coordinates), ]
  # subset species  from  our data
    DB_data <- read_excel("C:/Users/Dell-PC/Dropbox/CO_DBdata/abundance/Scarabaeinae_database_2024.xlsx", sheet="Scarabaeinae_database_2024")
    sp2 <- data.frame(subset(DB_data, scientificName == "Canthon_fulgidus_martinezi"))
    coordinates1 <- sp2[,c("lon_all_points", "lat_all_points")]
    coordinates1 <- na.omit(coordinates1)
    coordinates1 <- coordinates1[!duplicated(coordinates1), ]
    coordinates1$lon_all_points <- as.numeric(coordinates1$lon_all_points) 
    coordinates1$lat_all_points <- as.numeric(coordinates1$lat_all_points) 
#
  ggplot() +
    geom_sf(data = south_america, fill = "lightgrey", color = "black") +
    geom_point(data = coordinates, aes(x = decimalLongitude, y = decimalLatitude), color = "darkgreen", size = 3) +
    geom_point(data = coordinates1, aes(x = lon_all_points, y = lat_all_points), color = "darkblue", size = 3)
    