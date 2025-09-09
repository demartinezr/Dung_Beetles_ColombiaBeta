library(readxl)
library(mFD)
library(dplyr)
library(tidyr)
library(tibble)

# Overview of the functional framework

# load traits data
beha <- read_excel("C:/Users/PC/Dropbox/CO_DBdata/traits/behaviour/DB_Distributions_traits_2024.xlsx", sheet="DB_Distributions_traits")
beha <- beha[c("scientificName", "nest_guild", "diet_range", "activity")]
morpho <- read.table("C:/Users/PC/Dropbox/CO_DBdata/traits/morphometrics/morphometrics_mean.txt", header = TRUE)
morpho <- morpho[morpho$scientificName %in% beha$scientificName, ]

traits_DB <- as.data.frame(left_join(beha, morpho, by = "scientificName"))
rownames(traits_DB) <- traits_DB$scientificName
traits_DB <- traits_DB %>% select(-scientificName)

# load abundance data
db5 <- readRDS("C:/Users/PC/Dropbox/CO_DBdata/abundance/db5_distance.RDS")
db5$elev_standard_squared <- db5$elev_standard^2
db5$subregion_species <- paste0(db5$subregion, "__", db5$scientificName)
db5$cluster_species <- paste0(db5$cluster, "__", db5$scientificName)
db5$distance_from_range_scaled2 <- as.vector(db5$distance_from_range_scaled2)
# sum abundance day by species/point. Remove day covariate -> 241948 obs, 958 points, 243 species
db6 <- db5 %>%  
  group_by(across(-c(day, p_d))) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")  
# remove palm and young forest 219452 obs, 870 points, 243 species
db6 <- db6[db6$habitat !="Sy",]
db6 <- db6[db6$habitat !="PALM",]
db7 <- db6 |>
  filter(combinations ==1)

db7 <- db7 %>% filter(point %in% c("ALF1", "ALF2", "ALF3", "ALP1", "ALP2", "ALP3"))
db7 <- db7 %>% select(scientificName, abundance, habitat)
db8 <- db7 %>%
  group_by(habitat, scientificName) %>%
  summarise(total_abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(
    names_from = scientificName,
    values_from = total_abundance,
    values_fill = 0
  ) %>%
  column_to_rownames(var = "habitat")

################################################################################

# 1. Know your data

# 1.1. What types of traits am I using?
  
traits_DB_SM <- traits_DB %>%
  filter(rownames(traits_DB) %in% unique(db7$scientificName)) %>%
  mutate(
    nest_guild = as.factor(nest_guild),
    diet_range = as.factor(diet_range),
    activity   = as.factor(activity)
  )

traits_cat <- data.frame(
  trait_name = c("nest_guild", "diet_range", "activity", "bodysize", "legratio"),
  trait_type = c("N", "N", "N", "Q", "Q"),
  stringsAsFactors = FALSE
)

## 1.2. Summarize my traits

# Species traits summary:
traits_summ <- mFD::sp.tr.summary(
  tr_cat     = traits_cat,   
  sp_tr      = traits_DB_SM, 
  stop_if_NA = TRUE)

traits_summ$"tr_types" # Traits types 
traits_summ$"mod_list" # Traits types for non-continuous and non-fuzzy traits

# 1.3. Summarize my assemblages

# Summary of the assemblages * species dataframe:
db8 <- as.matrix(db8)
asb_sp_summ <- mFD::asb.sp.summary(asb_sp_w = db8)
head(asb_sp_summ$"asb_sp_occ", 3) 


asb_sp_occ <- asb_sp_summ$"asb_sp_occ"
asb_sp_summ$"sp_tot_w" # Species total biomass in all assemblage
asb_sp_summ$"asb_tot_w" # Total biomass per assemblage
asb_sp_summ$"asb_sp_richn" # Species richness per assemblage
asb_sp_summ$"asb_sp_nm"[[2]] # Names of species present in the first assemblage

#
# Gathering species into functional entities
# Computing distances between species based on functional traits

sp_dist <- mFD::funct.dist(
  sp_tr         = traits_DB_SM,
  tr_cat        = traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)

round(sp_dist, 3) 

# 4. Computing functional spaces & their quality
# 4.1. Compute multimensional functional spaces and assess their quality
fspaces_quality <- mFD::quality.fspaces(
  sp_dist             = sp_dist,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

round(fspaces_quality$"quality_fspaces", 3) # Quality metrics of spaces

# 4.2. Illustrating the quality of the selected functional spaces
mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality,
  quality_metric             = "mad",
  fspaces_plot               = c("tree_average", "pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

# 5. Test correlation between functional axes and traits
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

tr_faxes <- mFD::traits.faxes.cor(
  sp_tr          = traits_DB_SM, 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")], 
  plot           = TRUE)

# Print traits with significant effect:
tr_faxes$"tr_faxes_stat"[which(tr_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]
# print plot
tr_faxes$"tr_faxes_plot"

# 6. Plot functional space
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  color_bg        = "grey95",
  color_pool      = "darkgreen",
  fill_pool       = "white",
  shape_pool      = 21,
  size_pool       = 1,
  plot_ch         = TRUE,
  color_ch        = "black",
  fill_ch         = "white",
  alpha_ch        = 0.5,
  plot_vertices   = TRUE,
  color_vert      = "blueviolet",
  fill_vert       = "blueviolet",
  shape_vert      = 23,
  size_vert       = 1,
  plot_sp_nm      = NULL,
  nm_size         = 3,
  nm_color        = "black",
  nm_fontface     = "plain",
  check_input     = TRUE)

print(big_plot)


# Here are the plots for the fruits & baskets dataset for the first ten PCoA axis:
  
  big_plot <- mFD::funct.space.plot(
    sp_faxes_coord  = sp_faxes_coord,
    faxes           = NULL,
    name_file       = NULL,
    faxes_nm        = NULL,
    range_faxes     = c(NA, NA),
    color_bg        = "grey95",
    color_pool      = "darkgreen",
    fill_pool       = "white",
    shape_pool      = 21,
    size_pool       = 1,
    plot_ch         = TRUE,
    color_ch        = "black",
    fill_ch         = "white",
    alpha_ch        = 0.5,
    plot_vertices   = TRUE,
    color_vert      = "blueviolet",
    fill_vert       = "blueviolet",
    shape_vert      = 23,
    size_vert       = 1,
    plot_sp_nm      = NULL,
    nm_size         = 3,
    nm_color        = "black",
    nm_fontface     = "plain",
    check_input     = TRUE)

# Plot the graph with all pairs of axes:
big_plot$patchwork

# 7. Compute functional diversity indices & plot them
# 7.1. Functional alpha diversity indices in a multidimensional space

alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = db8,
  ind_vect         = c("fdis", "fmpd", "fnnd", "feve", "fric", "fdiv", "fori", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

fd_ind_values <- alpha_fd_indices$"functional_diversity_indices"
fd_ind_values

details_list <- alpha_fd_indices$"details"

plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("Forest", "Pasture"),
  ind_nm                   = c("fdis", "fide", "fnnd", "feve", "fric", 
                               "fdiv", "fori", "fspe"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_bg                 = "grey95",
  shape_sp                 = c(pool = 3, asb1 = 21, asb2 = 21),
  size_sp                  = c(pool = 0.7, asb1 = 1, asb2 = 1),
  color_sp                 = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_vert               = c(pool = "grey50", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_sp                  = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_vert                = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  color_ch                 = c(pool = NA, asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  fill_ch                  = c(pool = "white", asb1 = "#1F968BFF", asb2 = "#DCE319FF"),
  alpha_ch                 = c(pool = 1, asb1 = 0.3, asb2 = 0.3),
  shape_centroid_fdis      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fdiv      = c(asb1 = 22,  asb2 = 24),
  shape_centroid_fspe      = 23,
  color_centroid_fspe      = "black",
  size_sp_nm               = 3, 
  color_sp_nm              = "black",
  plot_sp_nm               = NULL,
  fontface_sp_nm           = "plain",
  save_file                = FALSE,
  check_input              = TRUE) 

plots_alpha$"fric"$"patchwork"
plots_alpha$"fdiv"$"patchwork"
plots_alpha$"fspe"$"patchwork"
plots_alpha$"fdis"$"patchwork"
plots_alpha$"fide"$"patchwork"
plots_alpha$"feve"$"patchwork"