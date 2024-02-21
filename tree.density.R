#### TREE DENSITY INDEX ####

#two subfunctions in Tree Density

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


# using *_ground_mask and *_canopy_mask to calculate densities
# first define function for using in moving window
tree_density <- function(dem, dsm, layer.code, size, output.dir){
  chm <- dsm-dem #canopy height model
  writeRaster(chm, paste0(output.dir,layer.code,"_chm.tif"),format="GTiff", overwrite=TRUE)
  veg <- chm>1 #above 1meter
  writeRaster(veg, paste0(output.dir,layer.code,"_veg.tif"),format="GTiff", overwrite=TRUE)
  # calculate total non zero pixels
  
  # density calculation 
  windowSize <- matrix(1, nrow = size, ncol = size) #matrix with identical weight and variable size of window
  # since the calculation is done on a binary raster, the result of a mean meanFunction is identical with canopyRatio function 
  # mean_chm <- focal(veg, w = windowSize, fun = meanFunction)
  # plot(chm_density)
  # writeRaster(mean_chm, paste0(output.dir, layer.code,"_",size, "_mean_above1m.tif"),format="GTiff", overwrite=TRUE)
  
  chm_density <- focal(veg, w = windowSize, fun = canopyRatio)
  #plot(chm_density)
  writeRaster(chm_density, paste0(output.dir, layer.code, "_",size,"_chm_density.tif"),format="GTiff", overwrite=TRUE)
  outputs <- list(
    chm = chm,         # canopy height
    veg = veg,       # canopy higher than 1m
    # mean_chm = mean_chm, # mean number of pixels above 1m in radius of moving window
    chm_density = chm_density# canopy density in radius in moving window
  )
  return(outputs)
}