library(raster)
library(RSAGA)
library(sf)
library(terra)
setwd("D:/Drone_DATA/2023_Oppkuven/")
env <- rsaga.env(workspace = "D:/Drone_DATA/2023_Oppkuven/",
path = "C:/Program Files/saga-7.4.0_x64") 
rsaga.get.libraries(path = env$modules)
rsaga.get.modules(libs = "ta_lighting")
####
# 10cm resolution
#### MASKING TO EXTENT ####
# example data 
#x <- lapply(ff, raster)
#t(sapply(x, function(i) as.vector(extent(i))))
# read GROUND in as text string
e1 <- "Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_e1/e1_02-ground_fill_0.1.tif"
e2 <- "Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_e2/e2_02-ground_fill_0.1.tif"
w1 <- "Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_w1/w1_02-ground_fill_0.1.tif"
w2 <- "Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_w2/w2_02-ground_fill_0.1.tif"
#
e1_g <- raster(e1)
e2_g <- raster(e2)
w1_g <- raster(w1)
w2_g <- raster(w2)

# MASK GROUND to minimize redundancy 
oppkuven_mask_e1 <- st_read("D:/Spatial_DATA/GEco - Projects/2023_Oppkuven_Msc_Ina_Mari/Oppkuven_QFIELD/oppkuven_area.gpkg", layer = "e1")
oppkuven_mask_e2 <- st_read("D:/Spatial_DATA/GEco - Projects/2023_Oppkuven_Msc_Ina_Mari/Oppkuven_QFIELD/oppkuven_area.gpkg", layer = "e2")
oppkuven_mask_w1 <- st_read("D:/Spatial_DATA/GEco - Projects/2023_Oppkuven_Msc_Ina_Mari/Oppkuven_QFIELD/oppkuven_area.gpkg", layer = "w1")
oppkuven_mask_w2 <- st_read("D:/Spatial_DATA/GEco - Projects/2023_Oppkuven_Msc_Ina_Mari/Oppkuven_QFIELD/oppkuven_area.gpkg", layer = "w2")
oppkuven_mask_e1<- st_zm(oppkuven_mask_e1, drop=TRUE, what = "ZM")# drop dimension Z from feature dataset
e1_ground_mask <- mask(e1_g, oppkuven_mask_e1)
e2_ground_mask <- mask(e2_g, oppkuven_mask_e2)
w1_ground_mask <- mask(w1_g, oppkuven_mask_w1)
w2_ground_mask <- mask(w2_g, oppkuven_mask_w2)

e1_ground_mask <- crop(e1_ground_mask, oppkuven_mask_e1)
e2_ground_mask <- crop(e2_ground_mask, oppkuven_mask_e2)
w1_ground_mask <- crop(w1_ground_mask, oppkuven_mask_w1)
w2_ground_mask <- crop(w2_ground_mask, oppkuven_mask_w2)

#it could be beneficial to implement CROP as well
writeRaster(e1_ground_mask, "01_masked/e1_ground_0.1.tif")
writeRaster(e2_ground_mask, "01_masked/e2_ground_0.1.tif")
writeRaster(w1_ground_mask, "01_masked/w1_ground_0.1.tif")
writeRaster(w2_ground_mask, "01_masked/w2_ground_0.1.tif", overwrite=TRUE)

# read CANOPY in as text string
e1_canopy <- "Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_e1/e1_03-canopy_0.1.tif"
e2_canopy <- "Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_e2/e2_03-canopy_0.1.tif"
w1_canopy <- "Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_w1/w1_03-canopy_0.1.tif"
w2_canopy <- "Oppkuven_PROCESSED_lidar/2023_26_10_Oppkuven_w2/w2_03-canopy_0.1.tif"

# read in as raster
e1_c <- raster(e1_canopy)
e2_c <- raster(e2_canopy)
w1_c <- raster(w1_canopy)
w2_c <- raster(w2_canopy)

# MASK CANOPY 
e1_canopy_mask <- mask(e1_c, oppkuven_mask_e1)
e2_canopy_mask <- mask(e2_c, oppkuven_mask_e2)
w1_canopy_mask <- mask(w1_c, oppkuven_mask_w1)
w2_canopy_mask <- mask(w2_c, oppkuven_mask_w2)
# crop to 
e1_canopy_mask <- crop(e1_canopy_mask, e1_ground_mask)
e2_canopy_mask <- crop(e2_canopy_mask, e2_ground_mask)
w1_canopy_mask <- crop(w1_canopy_mask, w1_ground_mask)
w2_canopy_mask <- crop(w2_canopy_mask, w2_ground_mask)

writeRaster(e1_canopy_mask, "01_masked/e1_canopy_0.1.tif")
writeRaster(e2_canopy_mask, "01_masked/e2_canopy_0.1.tif")
writeRaster(w1_canopy_mask, "01_masked/w1_canopy_0.1.tif")
writeRaster(w2_canopy_mask, "01_masked/w2_canopy_0.1.tif")


# check extent
extent(w2_ground_mask)
extent(w2_canopy_mask)
compareRaster(w2_ground_mask,w2_canopy_mask)
# par(mfrow=c(1,3))
# plot(w2_g)
# plot(oppkuven_mask_w2, add=TRUE)
# plot(w2_ground_crop)
# plot(oppkuven_mask_w2, add=TRUE)
# plot(w2_ground_mask)
# plot(oppkuven_mask_w2, add=TRUE)
stack(w2_ground_mask, w2_canopy_mask)
# check extent where x is the list of rasters
ff <- c(w2_ground_mask, w2_canopy_mask)
ff <- c(e1_ground_mask, e1_canopy_mask,
        e2_ground_mask, e2_canopy_mask,
        w1_ground_mask, w1_canopy_mask,
        w2_ground_mask, w2_canopy_mask)
list_of_rasters <- lapply(ff, raster) # where ff is array of strings
t(sapply(list_of_rasters, function(i) as.vector(extent(i))))# check extent
t(sapply(list_of_rasters, function(i) as.vector(res(i))))# check resolution
t(sapply(list_of_rasters, function(i) as.vector(origin(i))))# check origin

# resampling to keep extent https://stackoverflow.com/questions/65089257/how-to-keep-same-extent
rsts <- list(w2_ground_mask, w2_canopy_mask) #list(r1, r2, r3, r4)
for (i in 2:length(rsts)) {
  rsts[[i]] <- resample(rsts[[i]], r1)
} 
s <- stack(rsts)
s
# specific alignment 
# e2 and w1 still had different extent and origin
  

# ideally wrap into function
align_raster <- function(raster.file1, raster.file2){
  
}

c(ncol(e2_ground_mask),nrow(e2_ground_mask))
c(ncol(e2_canopy_mask),nrow(e2_canopy_mask))
origin(e2_ground_mask)
origin(e2_canopy_mask)
extent(e2_ground_mask)
extent(e2_canopy_mask)
# Calculate shift required
x_shift <- origin(e2_ground_mask)[1] - origin(e2_canopy_mask)[1]
# Shift the e2_canopy_mask raster
e2_canopy_mask <- shift(e2_canopy_mask, dx=x_shift, dy=0)
e2_extent <- intersect(extent(e2_ground_mask), extent(e2_canopy_mask))
e2_extent<- extent(584948,  585135, 6661632 ,6661860 )
e2_ground_mask <- crop(e2_ground_mask, e2_extent)
e2_canopy_mask <- crop(e2_canopy_mask, e2_extent)
# aligned_extent <- alignExtent(e2_canopy_mask, e2_ground_mask)
# raster2_aligned <- crop(e2_canopy_mask, aligned_extent)
# e2_ground_mask <- crop(e2_ground_mask, aligned_extent)

compareRaster(e2_ground_mask, e2_canopy_mask)


# e2 and w1 still had different extent origin
c(ncol(w1_ground_mask),nrow(w1_ground_mask))
c(ncol(w1_canopy_mask),nrow(w1_canopy_mask))
origin(w1_ground_mask)
origin(w1_canopy_mask)
extent(w1_ground_mask)
extent(w1_canopy_mask)
# Calculate shift required
x_shift <- origin(w1_ground_mask)[1] - origin(w1_canopy_mask)[1]
# Shift the e2_canopy_mask raster
w1_canopy_mask <- shift(w1_canopy_mask, dx=x_shift, dy=0)

w1_extent <- intersect(extent(w1_ground_mask), extent(w1_canopy_mask))
w1_extent<- extent(583468,  583579 , 6662088  ,6662309  )
w1_ground_mask <- crop(w1_ground_mask, w1_extent)
w1_canopy_mask <- crop(w1_canopy_mask, w1_extent)
compareRaster(w1_canopy_mask, w1_ground_mask)


# w2 still had different origin
c(ncol(w2_ground_mask),nrow(w2_ground_mask))
c(ncol(w2_canopy_mask),nrow(w2_canopy_mask))
origin(w2_ground_mask)
origin(w2_canopy_mask)
extent(w2_ground_mask)
extent(w2_canopy_mask)
compareRaster(w2_canopy_mask, w2_ground_mask)
# Calculate shift required
x_shift <- origin(w2_ground_mask)[1] - origin(w2_canopy_mask)[1]
# Shift the w2_canopy_mask raster
w2_canopy_mask <- shift(w2_canopy_mask, dx=x_shift, dy=0)
compareRaster(w2_canopy_mask, w2_ground_mask)

# w2_extent <- intersect(extent(w2_ground_mask), extent(w2_canopy_mask))
# w2_extent<- extent(584125.6 ,  585135, 6661632 ,6661860 )
# w2_ground_mask <- crop(w2_ground_mask, w2_extent)
# w2_canopy_mask <- crop(w2_canopy_mask, w2_extent)
# aligned_extent <- alignExtent(e2_canopy_mask, e2_ground_mask)
# raster2_aligned <- crop(e2_canopy_mask, aligned_extent)
# e2_ground_mask <- crop(e2_ground_mask, aligned_extent)

# e1_canopy_mask <- resample(e1_ground_mask, e1_canopy_mask)
# e2_canopy_mask <- resample(e2_ground_mask, e2_canopy_mask)
# w1_canopy_mask <- resample(w1_ground_mask, w1_canopy_mask)
# w2_canopy_mask <- resample(w2_ground_mask, w2_canopy_mask)


#### RESAMPLE THE ELEVATION TO 1, 5, 9M RESOLUTION

resample_raster <- function(tiff_path, resolution) {
  # Load the original raster
  original_raster <- raster(tiff_path)
  
  # Calculate the resampling factor
  original_res <- res(original_raster)[1] # assuming square pixels
  resample_factor <- resolution / original_res
  
  # Resample the raster
  resampled_raster <- aggregate(original_raster, fact=resample_factor)
  
  # write raster
  print(tiff_path)
  layer.name <- strsplit(tiff_path, split = "/")[[1]][3]
  layer.code <- paste0(strsplit(layer.name, split = "_")[[1]][1]) #split string
  writeRaster(resampled_raster, paste0("02_resampled/",layer.code,"_dem_", resolution,"m.tif" ),format="GTiff")
  # Return the resampled raster
  print(paste0("02_resampled/",layer.code,"_", resolution,".tif" ))
  return(resampled_raster)
  
}
# run the resampling for three resolutions ####
# resampled_raster <- resample_raster(w1, 9) # test

res <- c(1,5,9)
for (x in res){
  resample_raster(e1, x)
  }
for (x in res){
  resample_raster(e2, x)
}
for (x in res){
  resample_raster(w1, x)
}
for (x in res){
  resample_raster(w2, x)
}


#### CREATE DERIVED DATA FROM ELEVATION MODEL ####

# read all raster files inside HOME folder and add them to a list
# list of file names - strings
# areas <- list(e1, e2, w1, w2)
areas <- list.files(paste0(getwd(),'/01_masked/'), pattern="ground_0.1.tif", full.names=TRUE)
# read in raster layers and stack them into raster stack
#r.stack.cont <- stack(r.list.cont)
# basename(r.list.cont)

# structure for the loop of the FOUR areas
# file.name <- areas[4]
for (file.name in areas){
  print(file.name)
  layer.name <- strsplit(file.name, split = "/")[[1]][5]
  layer.code <- paste0(strsplit(layer.name, split = "_")[[1]][1])

      
      # derive topo variables in RSAGA
      rsaga.slope.asp.curv(in.dem = file.name,
                                    out.slope  = paste0("03_derived/", layer.code, "_slope.tif"),
                                    out.aspect = paste0("03_derived/", layer.code, "_aspect.tif"),
                                    out.cgene =  paste0("03_derived/", layer.code, "_cgene.tif"),
                                    out.cprof =  paste0("03_derived/", layer.code, "_cprof.tif"),
                                    out.cplan =  paste0("03_derived/", layer.code, "_cplan.tif"),
                                    unit.slope = "degrees",
                                    unit.aspect = "degrees", env = env)
      #https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.slope.asp.curv
      
      
      rsaga.wetness.index(in.dem = file.name,
                          out.wetness.index = paste0("03_derived/", layer.code, "_swi.tif"),
                          env = env)
      # Details The SAGA Wetness Index is similar to the Topographic Wetness Index (TWI), but it is based on a modified catchment area calculation (out.mod.carea), which does not treat the flow as a thin film as done in the calculation of catchment areas in conventional algorithms. As a result, the SWI tends to assign a more realistic, higher potential soil wetness than the TWI to grid cells situated in valley floors with a small vertical distance to a channel.
      # https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.wetness.index
      
      a.hillshade <- rsaga.hillshade(in.dem = file.name,
                                     out.grid = paste0("03_derived/", layer.code, "_hillshade.tif"),
                                     method = "standard",
                                     azimuth = 315,
                                     declination = 45,
                                     exaggeration = 4,
                                     env = env)
      #The Analytical Hillshading algorithm is based on the angle between the surface and the incoming light beams, measured in radians.
      #https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.hillshade
      
      # rsaga.insolation(in.dem = area.name,
      #                  out.direct = paste0("SAGA/", layer.code, "_insolation.direct.tif"),
      #                  out.diffuse = paste0("SAGA/", layer.code, "_insolation.diffuse.tif"),
      #                  out.total = paste0("SAGA/", layer.code, "_insolation.total.tif"),
      #                  time.step = 4, day.step = 1, days = 172:264, 
      #                  env = env)
      
      #in the Northern Hemisphere, summer might be considered from around June 21st (day 172) to September 21st (day 264).
      #https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.insolation
      # Calculation of incoming solar radiation (insolation). Based on the SADO (System for the Analysis of Discrete Surfaces) routines developed by Boehner & Trachinow.
      
      rsaga.pisr2(in.dem = file.name,
                  out.direct.grid = paste0("03_derived/", layer.code, "_insolation.direct.tif"),
                  out.diffuse.grid = paste0("03_derived/", layer.code, "_insolation.diffuse.tif"),
                  out.total.grid = paste0("03_derived/", layer.code, "_insolation.total.tif"),
                  out.duration =  paste0("03_derived/", layer.code, "_insolation.duration.tif"),
                  start.date = list(day=21,month=7,year=2023), 
                  end.date = list(day=21,month=10,year=2023),
                  time.step = 6,
                  day.step = 7,
                  env = env
                  )
      # 
      # https://rdrr.io/github/r-spatial/RSAGA/man/rsaga.pisr2.html
      # alternative paramters from https://gis.stackexchange.com/questions/307426/rsaga-error-using-rsaga-pisr2
      # location = "grid", latitude = 48.5, unit = "kWh/m2", solconst = 1367, method = "height", hgt.atmosphere = 12000, cmp.pressure = 1013, cmp.water.content = 1.68, cmp.dust = 100, lmp.transmittance = 70, time.range = c(0, 24), time.step = 0.5, start.date = list(day = 1, month = 2, year = 2015), end.date = list(day = 2, month = 2, year = 2015), day.step = 5,
      # catchment area
      rsaga.topdown.processing(in.dem = file.name,
                               out.carea = paste0("03_derived/", layer.code, "_catch_area.tif"),
                               method = "mfd", 
                               env = env)
      
      #https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.topdown.processing
      #Calculate the size of the local catchment area (contributing area), accumulated material, and flow path length, using top-down processing algorithms from the highest to the lowest cell.
}



#### TREE DENSITY INDEX ####


# using *_ground_mask and *_canopy_mask to calculate densities
# first define function for using in moving window
meanFunction <- function(x) { mean(x, na.rm = TRUE) } 
# 
canopyRatio <- function(matrix) {
  canopyCells <- sum(matrix == 1, na.rm = TRUE)  # Count canopy cells
  nonCanopyCells <- sum(matrix == 0, na.rm = TRUE)  # Count non-canopy cells
  
  if (nonCanopyCells == 0) {
    return(1)  # All cells are canopy cells
  } else {
    return(canopyCells / (canopyCells + nonCanopyCells))  # Ratio
  }
}
## tree density 
locality <- c("e1", "e2", "w1", "w2")
#locality[1]

# define the three spatial resolutions 
a=1/res(e1_ground_mask) # how many pixels(vertical,horizontal) to take into focal account to achieve 1m resolution, 3m and 9m,
b=5/res(e1_ground_mask)
c=9/res(e1_ground_mask)

size <- list(a,b,c)
size[[1]][1] # first element from the first list

# windowSize <- matrix(1, nrow = 9, ncol = 9)
# windowSize <- matrix(1, nrow = 49, ncol = 49)
# windowSize <- matrix(1, nrow = 89, ncol = 89)

#steps in tree density calculations

tree_density <- function(dem, dsm, locality, size){
    chm <- dsm-dem #canopy height model
    writeRaster(chm, paste0("03_derived/",locality,"_chm.tif"),format="GTiff", overwrite=TRUE)
    veg_1 <- chm>1 #above 1meter
    writeRaster(veg_1, paste0("03_derived/",locality,"_veg_above1m.tif"),format="GTiff", overwrite=TRUE)
    # calculate total non zero pixels
    
    # density calculation 
    windowSize <- matrix(1, nrow = size, ncol = size) #matrix with identical weight and variable size of window
    mean_chm <- focal(veg_1, w = windowSize, fun = meanFunction)
    #plot(chm_density)
    writeRaster(mean_chm, paste0("03_derived/", locality,"_",size, "_mean_above1m.tif"),format="GTiff", overwrite=TRUE)
    
    chm_density_1 <- focal(veg_1, w = windowSize, fun = canopyRatio)
    #plot(chm_density_1)
    writeRaster(chm_density_1, paste0("03_derived/", locality, "_",size,"_chm_density_above1m.tif"),format="GTiff", overwrite=TRUE)
    outputs <- list(
      output1 = chm,         # canopy height
      output2 = veg_1,       # canopy higher than 1m
      output3 = mean_chm, # mean number of pixels above 1m in radius of moving window
      output4 = chm_density_1# canopy density in radius in moving window
      )
    return(outputs)
    }
# run and save into variable with four results
# canopy_analysis <- tree_density(dem=testing1[[5]], dsm = testing2[[5]], locality = locality[1], size = 9)
# canopy_analysis <- tree_density(dem=w2_ground_mask, dsm = w2_canopy_mask, locality = "w2" , size = 9)

#locality <- c("e1", "e2", "w1", "w2")
locality <- c( "w2")
for (x in locality){
  tree_density(dem=eval(parse(text=paste0(x,"_ground_mask"))), 
               dsm=eval(parse(text=paste0(x,"_canopy_mask"))), 
               locality = x, 
               size = 9)
  
  tree_density(dem=eval(parse(text=paste0(x,"_ground_mask"))), 
               dsm=eval(parse(text=paste0(x,"_canopy_mask"))),
               locality = x, 
               size = 49)
  
  tree_density(dem=eval(parse(text=paste0(x,"_ground_mask"))), 
               dsm=eval(parse(text=paste0(x,"_canopy_mask"))),
               locality = x, 
               size = 89)
}
# do this for all sites 

# 

# 
par(mfrow=c(1,4))
    plot(canopy_analysis$output1) 
    plot(canopy_analysis$output2) 
    plot(canopy_analysis$output3) 
    plot(canopy_analysis$output4) 
   
    
# possible to use loop for different resolutions

res <- c(1,5,9)
for (x in res){
  resample_raster(e1, x)
}
#### STACK ALL CREATED RASTERS ####
# check which raster has the smallest extent
# crop to this extent
# stack


#### read in student data ####

ina_pts <- st_read("Ina_sites.gpkg")
mari_pts <- st_read("Mari_sites.gpkg")

#### EXTRACT VALUES TO POINTS ####

# W2 refers to INA's data 
# E1/E2/W1/W2 to Mari's data 
# load in csv data
# pts <- read.csv2("C:/Users/peterhor/OneDrive - Universitetet i Oslo/2019 - Natural History Museum/MSc_PhD_Supervision/2022_Ina/80_precise_coordinates.csv")
# pts$X <-as.numeric(pts$X)
# pts$Y <-as.numeric(pts$Y)
# pts_utm <- SpatialPoints(cbind(pts$X,pts$Y))

# INA
w2_layer_list <- list.files(paste0(getwd(),'/03_derived/'), pattern="w2_", full.names=TRUE)
# read in raster layers and stack them into raster stack
w2_layers <- stack(w2_layer_list)

ina_pts_values <- extract(w2_layers, ina_pts)
ina_pts_t <- cbind(ina_pts, ina_pts_values)
write.csv2(ina_pts_t, paste0(getwd(),'/05_output/'))

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

