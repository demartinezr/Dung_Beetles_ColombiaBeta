# This script was devolped by Jacob socolar for birds and adapted for dung
# beetles by DEMR
library(raster)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
llanos_amazon <- sf::st_transform(sf::st_read("/Users/jacobsocolar/Dropbox/Work/Colombia/Data/GIS/masking_polygons/E_Amazon_Llanos.kml"),
                                  AEAstring)