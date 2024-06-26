library(raster)
library(RSAGA)
library(sf)
library(terra)
DataDir <- "D:/Drone_DATA/2023_Oppkuven/"

# call functions and execute


# list_of_rasters <- lapply(c(dtm_path,dsm_path),  raster) # where we load  array of strings
# t(sapply(list_of_rasters, function(i) as.vector(extent(i))))# check extent
# t(sapply(list_of_rasters, function(i) as.vector(res(i))))# check resolution
# t(sapply(list_of_rasters, function(i) as.vector(origin(i))))# check origin


### For each locality ####
locality <- c("e1", "e2", "w1", "w2")
for (layer.code in locality){
  message(layer.code)
  
  dtm_path <- paste0(DataDir, "Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_", layer.code, "/", layer.code, "_02-ground_fill_0.1.tif")
  dsm_path <- paste0(DataDir, "Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_", layer.code, "/", layer.code, "_03-canopy_0.1.tif")
  #
  dtm <- raster(dtm_path)
  dsm <- raster(dsm_path)
  
  #### MAKE STACKABLE ####
  # import funciton from file 
  source("make.stackable.R")
  
  
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
