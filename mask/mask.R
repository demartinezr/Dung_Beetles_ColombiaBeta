# This script generates a mask for the regions of colombia that we do not analyze.

library(sf)
library(raster)
library(ggplot2)

`%ni%` <- Negate(`%in%`)
AEAstring <- "+proj=aea +lat_1=-4.2 +lat_2=12.5 +lat_0=4.1 +lon_0=-73 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

source("D:/Doctorado/Tesis/Topographic_Col/unidades_topograficas.R")
pacific <- st_transform(pacific, crs(raster_elev_AEA))
tacarcuna <- st_transform(tacarcuna, crs(raster_elev_AEA))
snsm <- st_transform(snsm, crs(raster_elev_AEA))

raster_elev_AEA <- raster("D:/Capas/America/dem/elev_raster/raster_elev_AEA.grd")

pac_raster <- fasterize::fasterize(pacific, raster_elev_AEA)
pac_raster[raster_elev_AEA > 1100] <- NA
mask_raster <- pac_raster
tac_raster <- fasterize::fasterize(tacarcuna, raster_elev_AEA)
mask_raster[tac_raster == 1] <- 1
snsm_raster <- fasterize::fasterize(snsm, raster_elev_AEA)
snsm_raster[raster_elev_AEA < 3000] <- NA
snsm_raster[is.na(raster_elev_AEA)] <- NA
mask_raster[snsm_raster == 1] <- 1
if (!requireNamespace("colmaps", quietly = TRUE)) {
  devtools::install_github("nebulae-co/colmaps")
}
dpts <- colmaps::departamentos

mask_dpts <- dpts[dpts$depto %in% c("Cesar", "La Guajira", "Magdalena", "Bolívar", "Atlántico", "Córdoba", "Sucre", "Antioquia", "Chocó"), ]

md <- st_buffer(st_as_sf(mask_dpts), 3000)
md <- st_cast(md, "MULTIPOLYGON")

md_raster <- fasterize::fasterize(md, raster_elev_AEA)
md_raster[raster_elev_AEA > 1100] <- NA

mask_raster[md_raster == 1] <- 1
writeRaster(mask_raster, "D:/Capas/America/dem/mask/mask.grd", overwrite=T)
par(mar=c(2,2,1,1))
plot(mask_raster)
