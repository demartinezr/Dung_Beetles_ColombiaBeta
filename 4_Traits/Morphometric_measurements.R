# Protocol to unify morphometric measurements of 243 species/morpho species from 
# all projects, new pictures from 2022
#
setwd("C:/Users/Dell-PC/Dropbox/CO_DBdata")
#
# packages
library(readxl)
library(dplyr)
#
# dataset new photos morphometrics
  data_new <- read_excel("./traits/morphometrics/template_sheet_morph_traits_2024.xlsx", sheet="Sheet1")
  data_new$backtibialength <- as.numeric(data_new$backtibialength)
  data_new$backspurlength <- as.numeric(data_new$backspurlength)
# filter and mean by individual_id_number para elytrawidth, bodywidth, elytralength, bodylength
  mean_body_elytra <- data_new  %>%
  filter(grepl("_1|_2$", photo_number)) %>%
  group_by(scientificName, individual_id_number) %>%
  summarize(across(c(elytrawidth, bodywidth, elytralength, bodylength), mean, na.rm = TRUE), .groups = "drop")
# filter and mean by individual_id_number para frontfemurlength, fronttibialength, backfemurlength, backtibialength, backspurlength
  mean_legs <- data_new %>%
  filter(grepl("_3|_4$", photo_number)) %>%
  group_by(scientificName, individual_id_number) %>%
  summarize(across(c(frontfemurlength, fronttibialength, backfemurlength, backtibialength, backspurlength), mean, na.rm = TRUE), .groups = "drop")
# merge by individual_id_number y scientificName
  means_new <- merge(mean_body_elytra, mean_legs, by = c("scientificName", "individual_id_number"))
#
# dataset 2021 morphometrics
  data_2021 <- read.csv("./traits/morphometrics/output_DB_rawdata_alltraits_combinedMASTER.csv", header =TRUE, sep =";")
#
# Intersect colnames between data frames
  common_cols <- intersect(names(data_2021), names(means_new))
# select common columns between data frames
  data_2021_common <- data_2021 %>%
  select(all_of(common_cols))
  means_new_common <- means_new %>%
  select(all_of(common_cols))
# Join dataframes based on common columns
  morphometrics <- bind_rows(data_2021_common, means_new_common)
  morphometrics <- morphometrics %>%
  arrange(scientificName)
# traits body size and leg ratio measurements
  morphometrics$bodysize <- as.numeric(morphometrics$elytrawidth*
                                         morphometrics$bodylength)
  morphometrics$legratio <- as.numeric((morphometrics$frontfemurlength+
       morphometrics$fronttibialength)/
         (morphometrics$backfemurlength+
          morphometrics$backtibialength+
          morphometrics$backspurlength))
# mean id duplicates (messurements by DEMR and Felicity)  
  morphometrics <- morphometrics %>%
    group_by(individual_id_number) %>%
    summarise(
      scientificName = first(scientificName),
      across(c(elytrawidth, elytralength, bodylength, bodywidth,
                       frontfemurlength, fronttibialength, backfemurlength, 
                       backtibialength, backspurlength, bodysize, legratio), 
                     mean, na.rm = TRUE))
  write.table(morphometrics, "./traits/morphometrics/output_rawdata_alltraits_2024.txt")
  morphometrics <- read.table("./traits/morphometrics/output_rawdata_alltraits_2024.txt", header=TRUE)
  
  # traits body size and leg ratio means
  morphometrics_mean <- morphometrics %>%
    group_by(scientificName) %>%
    summarise(
      across(
        .cols = c(bodysize, legratio),
        .fns = ~ mean(.x[!is.na(.x) & !is.nan(.x) & !is.infinite(.x) & .x != 0], na.rm = TRUE),
        .names = "{col}"
      )
    ) %>%
    ungroup()
  write.table(morphometrics_mean, "./traits/morphometrics/morphometrics_mean.txt")
#
# Merge morphometrics measurements and behaviour traits in Dung beetles data set
# DB_data <- read_excel("C:./abundance/Scarabaeinae_database_2024.xlsx", sheet="Scarabaeinae_database_2024")
# behaviour <- read_excel("./traits/behaviour/DB_Distributions_traits_2024.xlsx", sheet="DB_Distributions_traits")
# behaviour <- behaviour[c("scientificName", "nest_guild", "diet_range", "activity")]
# DB_data <- left_join(DB_data, morphometrics_mean, by = "scientificName") 
# DB_data <- left_join(DB_data, behaviour, by = "scientificName")
 