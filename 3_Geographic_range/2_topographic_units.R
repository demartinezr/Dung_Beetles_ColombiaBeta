# ------------------------------------------------------------------------------
# This script defines topographic units of Colombia used for geographic range
# estimation. It combines HydroSHEDS basin polygons with manually defined
# biogeographic mountain boundaries to generate spatial units representing
# major drainage systems and mountain barriers.
#
# HydroSHEDS basin levels (2–5) are used to delineate large drainage regions
# (e.g., Amazon–Orinoco, Magdalena, Cauca, Pacific). These regions are then
# refined using mountain polygons (e.g., Andes ranges, Sierra Nevada de Santa
# Marta, Darién ranges) to subdivide basins according to major topographic
# barriers that influence species distributions.
#
# The resulting polygons represent topographic units used to constrain the
# estimation of species geographic ranges.
# ------------------------------------------------------------------------------
setwd("F:/Doctorado/Tesis/Topographic_Col/inputs")
library(ggplot2)
library(sf)
library(dplyr)

# Load HydroSHEDS basin layers at multiple hierarchical levels
# (levels 2–5) for South America
hs2 <- st_read(paste0('./hydrosheds/hybas_sa_lev01-06_v1c/hybas_sa_lev02_v1c.shp'))
hs3 <- st_read(paste0('./hydrosheds/hybas_sa_lev01-06_v1c/hybas_sa_lev03_v1c.shp'))
hs4 <- st_read(paste0('./hydrosheds/hybas_sa_lev01-06_v1c/hybas_sa_lev04_v1c.shp'))
hs5 <- st_read(paste0('./hydrosheds/hybas_sa_lev01-06_v1c/hybas_sa_lev05_v1c.shp'))
# Define major drainage regions by combining selected HydroSHEDS basins
amazon_orinoco <- st_make_valid(st_union(st_make_valid(hs2[2, ]), st_make_valid(hs3[4, ])))
pacific_prelim <- st_make_valid(st_union(st_make_valid(hs3[1,]), st_make_valid(hs3[25,])))
gen_magdalena <- st_make_valid(hs3[2, ])
cauca <- st_make_valid(hs5[5, ])
valledupar <- st_make_valid(hs5[8, ])
magdalena <- st_make_valid(st_difference(gen_magdalena, st_make_valid(st_union(cauca, valledupar))))[,c('HYBAS_ID', 'geometry')]
catatumbo <- st_make_valid(hs5[12, ])
maracaibo <- st_make_valid(hs3[3, ])
guajira_valledupar <- st_make_valid(st_union(valledupar, st_difference(maracaibo, catatumbo))[c('HYBAS_ID', 'geometry')])
# Load mountain polygons used to subdivide basins according to
# major topographic barriers (e.g., Andes, Sierra Nevada de Santa Marta)
central <- st_make_valid(st_zm(st_read(paste0('./biogeographic_clips/mountains_clips/centralAndes__cauca__magdalena.kml'))))
st_is_valid(central) # see https://github.com/r-spatial/sf/issues/1771
sf_use_s2(FALSE)
central <- st_make_valid(central)
sf_use_s2(TRUE)
st_is_valid(central)

snsm <- st_make_valid(st_zm(st_read(paste0('./biogeographic_clips/mountains_clips/SNSM__guajira_valledupar.kml'))))
pasto <- st_make_valid(st_zm(st_read(paste0('./biogeographic_clips/mountains_clips/pasto__pacific.kml'))))
tacarcuna <- st_make_valid(st_zm(st_read(paste0('./biogeographic_clips/mountains_clips/tacarcuna__pacific.kml'))))

# Refine basin boundaries by intersecting and subtracting mountain regions
# to generate topographically constrained spatial units
pacific_prelim2 <- st_make_valid(st_difference(pacific_prelim, pasto))
pacific <- st_make_valid(st_difference(pacific_prelim2, tacarcuna))
pasto <- st_make_valid(st_intersection(pacific_prelim, pasto))
tacarcuna <- st_make_valid(st_intersection(pacific_prelim, tacarcuna))

cauca_west <- st_make_valid(st_difference(cauca, central))
cauca_east <- st_make_valid(st_intersection(cauca, central))

magdalena_west <- st_make_valid(st_intersection(magdalena, central))
magdalena_east <- st_make_valid(st_difference(magdalena, central))

guajira_perija <- st_make_valid(st_difference(guajira_valledupar, snsm))
snsm <- st_make_valid(st_intersection(guajira_valledupar, snsm))

# Compile all resulting topographic units into a named list
out <- list(
  amazon_orinoco = amazon_orinoco,
  catatumbo = catatumbo,
  guajira_perija = guajira_perija,
  snsm = snsm,
  magdalena_east = magdalena_east,
  magdalena_west = magdalena_west,
  cauca_east = cauca_east,
  cauca_west = cauca_west,
  pacific = pacific,
  patia = pasto,
  tacarcuna = tacarcuna
)

# out_combined <- bind_rows(lapply(names(out), function(name) {
#  out[[name]] %>% mutate(region = name)
# }))

# divisions <- st_read("D:/Capas/Colombia/Colombia/COL_adm0.shp")
# topographic_cliped <- st_intersection(out_combined, divisions)

# saveRDS(topographic_cliped, "C:/Users/Dell-PC/Dropbox/CO_DBdata/SIG/topographic_units/topographic_units.rds")

# ggplot(data = topographic_cliped) +
#   geom_sf(aes(fill = region), color = "black") + 
#   scale_fill_viridis_d() +                       
#   theme_minimal() +                              
#   labs(title = "Topographic units", fill = "Region") +
#   theme_classic()
