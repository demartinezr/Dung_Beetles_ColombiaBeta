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
