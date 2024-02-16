
# resample_raster <- function(tiff_path, resolution) {
#   # Load the original raster
#   original_raster <- raster(tiff_path)
#   
#   # Calculate the resampling factor
#   original_res <- res(original_raster)[1] # assuming square pixels
#   resample_factor <- resolution / original_res
#   
#   # Resample the raster
#   resampled_raster <- aggregate(original_raster, fact=resample_factor)
#   
#   # write raster
#   print(tiff_path)
#   layer.name <- strsplit(tiff_path, split = "/")[[1]][3]
#   layer.code <- paste0(strsplit(layer.name, split = "_")[[1]][1]) #split string
#   writeRaster(resampled_raster, paste0("02_resampled/",layer.code,"_dem_", resolution,"m.tif" ),format="GTiff")
#   # Return the resampled raster
#   print(paste0("02_resampled/",layer.code,"_", resolution,".tif" ))
#   return(resampled_raster)
#   
# }

resample_raster_stack <- function(tiff_path, resolution, out_dir) {
  # Load the original raster stack
  original_stack <- stack(tiff_path)
  
  # Initialize an empty stack to store the resampled layers
  resampled_stack <- stack()
  
  # Iterate over each layer in the stack
  for (i in 1:nlayers(original_stack)) {
    # Select the current layer
    original_raster <- original_stack[[i]]
    
    # Calculate the resampling factor
    original_res <- res(original_raster)[1]  # Assuming square pixels
    resample_factor <- resolution / original_res
    
    # Resample the raster
    resampled_raster <- aggregate(original_raster, fact=resample_factor)
    
    # Add the resampled raster to the stack
    resampled_stack <- stack(resampled_stack, resampled_raster)
  }
  
  # Write rasters and prepare return information
  resampled_file_paths <- character(nlayers(resampled_stack))
  for (i in 1:nlayers(resampled_stack)) {
    layer.name <- names(resampled_stack)[i]
    output_path <- paste0(out_dir, resolution, "m_", layer.name, ".tif")
    writeRaster(resampled_stack[[i]], output_path, format="GTiff")
    resampled_file_paths[i] <- output_path
  }
  
  # Print information and return the resampled stack
  print(resampled_file_paths)
  return(resampled_stack)
}