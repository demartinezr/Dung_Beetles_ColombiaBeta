# ------------------------------------------------------------------------------
# This script extracts the mainland polygon of Colombia from the national
# administrative boundary shapefile (GADM). The original shapefile contains
# multiple disjoint polygons representing the mainland and offshore islands.
# To ensure consistent spatial analyses for geographic range estimation,
# the script identifies the largest polygon (mainland Colombia) and removes
# smaller island polygons. The resulting mainland boundary is used as the
# base spatial extent for subsequent geographic range calculations.
# ------------------------------------------------------------------------------

library(sf)
# read in GADM colombia shapefile
colombia <- st_read("F:/Capas/Colombia/Colombia/COL_adm0.shp")

# Calculate the number of vertices for each polygon to identify the largest one (mainland)
npoly <- length(colombia$geometry[[1]])
size <- rep(0,npoly)
for(i in 1:npoly){
  size[i] <- dim(colombia$geometry[[1]][[i]][[1]])[1]
}
# Select the largest polygon, corresponding to mainland Colombia
divisions <- colombia$geometry[[1]][[which(size == max(size))]] %>%
  st_polygon() %>% st_sfc() %>% st_sf()
st_crs(divisions) <- st_crs(divisions)

library(sf)
# read in GADM colombia shapefile
colombia <- st_read("F:/Capas/Colombia/Colombia/COL_adm0.shp")

# File consists of many disjoint polygons, representing the mainland and numerous islands. Figure out which is the mainland
# and extract it.
npoly <- length(colombia$geometry[[1]])
size <- rep(0,npoly)
for(i in 1:npoly){
  size[i] <- dim(colombia$geometry[[1]][[i]][[1]])[1]
}

divisions <- colombia$geometry[[1]][[which(size == max(size))]] %>%
  st_polygon() %>% st_sfc() %>% st_sf()
st_crs(divisions) <- st_crs(divisions)
