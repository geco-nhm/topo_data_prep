library(raster)
library(RSAGA)
library(sf)
library(terra)
setwd("D:/Drone_DATA/2023_Oppkuven/")

# call functions and execute


# list_of_rasters <- lapply(c(dtm_path,dsm_path),  raster) # where we load  array of strings
# t(sapply(list_of_rasters, function(i) as.vector(extent(i))))# check extent
# t(sapply(list_of_rasters, function(i) as.vector(res(i))))# check resolution
# t(sapply(list_of_rasters, function(i) as.vector(origin(i))))# check origin


### For each locality ####
locality <- c("e1", "e2", "w1", "w2")
for (layer.code in locality){
  
  dtm_path <- paste0("Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_", layer.code, "/", layer.code, "_02-ground_fill_0.1.tif")
  dsm_path <- paste0("Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_", layer.code, "/", layer.code, "_03-canopy_0.1.tif")
  #
  dtm <- raster(dtm_path)
  dsm <- raster(dsm_path)
  
  #### MAKE STACKABLE ####
  # import funciton from file 
  source("C:/Users/peterhor/Documents/GitHub/GIS_data_prep/make.stackable.R")
  
  
  r.stack <- make.stackable(dtm,dsm)
  
  plot(r.stack)
  
  
  # then proceed to masking by shapefile
  oppkuven_mask <- st_read("D:/Spatial_DATA/GEco - Projects/2023_Oppkuven_Msc_Ina_Mari/oppkuven_area.gpkg", layer = layer.code)
  oppkuven_mask<- st_zm(oppkuven_mask, drop=TRUE, what = "ZM")
  r.stack.crop <- crop(r.stack, oppkuven_mask)
  r.stack.mask <- mask(r.stack.crop, oppkuven_mask)
  # drop dimension Z from feature dataset
  
  plot(r.stack.mask)
  
  
  dtm_mask_path=paste0("01_masked/", layer.code, "_", names(r.stack.mask$dtm), ".tiff")
  writeRaster(r.stack.mask$dtm, filename = dtm_mask_path, format = "GTiff", overwrite=TRUE)
  dsm_mask_path=paste0("01_masked/", layer.code, "_", names(r.stack.mask$dsm), ".tiff")
  writeRaster(r.stack.mask$dsm, filename = dsm_mask_path, format = "GTiff", overwrite=TRUE)
  
  
  #### CREATE DERIVED DATA FROM ELEVATION MODEL ####
  source("C:/Users/peterhor/Documents/GitHub/GIS_data_prep/derive.vars.R")
  
  derive.vars(output.dir = "02_derived/", dtm = paste0(getwd(),"/", dtm_mask_path), layer.code = layer.code)
  
  
  #### TREE DENSITY INDEX ####
  source("C:/Users/peterhor/Documents/GitHub/GIS_data_prep/tree.density.R")
  # how many pixels(vertical,horizontal) to take into focal account to achieve 1m resolution, 3m and 9m,
  a=1/res(r.stack.mask$dtm) 
  b=5/res(r.stack.mask$dtm)
  c=9/res(r.stack.mask$dtm)
  
  size <- list(a,b,c)
  
  # run and save into variable with four results
  
  #dem, dsm, layer.code, size, output.dir
  tree_density(dem = r.stack.mask$dtm, dsm = r.stack.mask$dsm, layer.code = layer.code , size = 9, output.dir = "02_derived/")
  tree_density(dem = r.stack.mask$dtm, dsm = r.stack.mask$dsm, layer.code = layer.code , size = 49, output.dir = "02_derived/")
  tree_density(dem = r.stack.mask$dtm, dsm = r.stack.mask$dsm, layer.code = layer.code , size = 89, output.dir = "02_derived/")
  
  
  
  #### RESAMPLE ####
  
  # run the resampling for three resolutions 
  # resampled_raster <- resample_raster(w1, 9) # test
  source("C:/Users/peterhor/Documents/GitHub/GIS_data_prep/resample.raster.R")
  
  list.prepared <- list.files(paste0(getwd(),"/02_derived/"), pattern= paste0(layer.code, "_"), full.names=TRUE)
  
  resample_raster_stack(tiff_path = list.prepared, resolution = 1, out_dir="03_resampled/")
  resample_raster_stack(tiff_path = list.prepared, resolution = 5, out_dir="03_resampled/")
  resample_raster_stack(tiff_path = list.prepared, resolution = 9, out_dir="03_resampled/")
  

}
#### EXTRACT points ####
#### read in student data ####

ina_pts <- st_read("Ina_sites.gpkg")
mari_pts <- st_read("Mari_sites.gpkg")

# INA
# first take the original resolution 0.1m
# read in raster layers and stack them into raster stack
w2_layer_list <- list.files(paste0(getwd(),'/02_derived/'), pattern=paste0(layer.code, "_"), full.names=TRUE)
w2_layers <- stack(w2_layer_list)
ina_pts_values <- extract(w2_layers, ina_pts)
ina_pts_t <- cbind(ina_pts, ina_pts_values)

#then take the three resolutions at 1,5,9m
res <- c("1m","5m","9m")
for(x in res){
  w2_layer_list <- list.files(paste0(getwd(),'/03_resampled/'), pattern=paste0(x,"_",layer.code, "_"), full.names=TRUE)
  # read in raster layers and stack them into raster stack
  w2_layers <- stack(w2_layer_list)
  
  ina_pts_values <- extract(w2_layers, ina_pts)
  ina_pts_t <- cbind(ina_pts_t, ina_pts_values)
}
write.csv2(ina_pts_t, paste0(getwd(),'/05_output/ina_pts.csv'))


# MARI
#e1
e1_layer_list <- list.files(paste0(getwd(),'/03_derived/'), pattern="e1_", full.names=TRUE)
# read in raster layers and stack them into raster stack
e1_layers <- stack(e1_layer_list)
#e2
e2_layer_list <- list.files(paste0(getwd(),'/03_derived/'), pattern="e2_", full.names=TRUE)
# read in raster layers and stack them into raster stack
e2_layers <- stack(e2_layer_list)
#w1
w1_layer_list <- list.files(paste0(getwd(),'/03_derived/'), pattern="w1_", full.names=TRUE)
# read in raster layers and stack them into raster stack
w1_layers <- stack(w1_layer_list)
#w2
w2_layer_list <- list.files(paste0(getwd(),'/03_derived/'), pattern="w2_", full.names=TRUE)
# read in raster layers and stack them into raster stack
w2_layers <- stack(w2_layer_list)

mari_pts_values_e1 <- extract(e1_layers, mari_pts)
mari_pts_values_e2 <- extract(e2_layers, mari_pts)
mari_pts_values_w1 <- extract(w1_layers, mari_pts)
mari_pts_values_w2 <- extract(w2_layers, mari_pts)


mari_pts_t <- cbind(ina_pts, mari_pts_values_e1, mari_pts_values_e2, mari_pts_values_w1, mari_pts_values_w2)
write.csv2(mari_pts_t, paste0(getwd(),'/05_output/mari_pts.csv'))

