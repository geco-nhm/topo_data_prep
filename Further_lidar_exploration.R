library(lidR)
library(sf)
library(terra)
library(ggplot2)
setwd("D:/Drone_DATA/2023_Oppkuven")
las <- readLAS("Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_w1/w1_01-compressed.laz")
las@crs#"EPSG",32632      
print(las)

las_orig <- readLAS("Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_w1/w1_cloud.las")
Encoding(las_orig)
# read in polygon
oppkuven_mask <- st_read("D:/Spatial_DATA/GEco - Projects/2023_Oppkuven_Msc_Ina_Mari/oppkuven_area.gpkg", layer = "w1")
oppkuven_mask$geom
# it is in WGS 84 UTM 32, convert into ETRS 89 UTM 32
opp <- st_transform(oppkuven_mask, "EPSG:32632")

clipped_las = clip_roi(las_orig, opp)
print(clipped_las)
las_check(clipped_las)
# plot(las)

#read in mari points
mari_points <- st_read("Mari_sites.gpkg")
mari_points$geometry[c(1,5)]
p1 <- c(583509, 6662180)
p2 <- c(583536, 6662180)
las_tr <- clip_transect(clipped_las, p1, p2, width = 4, xz = TRUE)

ggplot(las_tr@data, aes(X,Z, color = Z)) + 
  geom_point(size = 0.5) + 
  coord_equal() + 
  theme_minimal() +
  scale_color_gradientn(colours = height.colors(50))


chm <- rasterize_canopy(clipped_las, 0.5, pitfree(subcircle = 0.2))
plot(clipped_las, bg = "white", size = 4)


ttops <- locate_trees(clipped_las, lmf(ws = 5))

plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

ttops_3m <- locate_trees(clipped_las, lmf(ws = 3))
ttops_11m <- locate_trees(clipped_las, lmf(ws = 11))

par(mfrow=c(1,2))
plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops_3m), add = TRUE, pch = 3)
plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops_11m), add = TRUE, pch = 3)