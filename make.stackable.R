

#TODO make for list of rasters
make.stackable <- function(dtm, dsm){
  # check Origin
  # Calculate shift required for X and Y origin coordinates 
  x_shift <- origin(dtm)[1] - origin(dsm)[1]
  y_shift <- origin(dtm)[2] - origin(dsm)[2]
  if (max(x_shift, y_shift) > 1){
    stop("Origin shift too high")
  }
    # Shift the surface raster
  dsm <- shift(dsm, dx=x_shift, dy=y_shift)
  print("canopy model shifted")
  
  
   # compare the extents 
  xmin(dtm)
  new_xmin = max(xmin(dtm), xmin(dsm))
  new_ymin = max(ymin(dtm), ymin(dsm))
  new_xmax = min(xmax(dtm), xmax(dsm))
  new_ymax = min(ymax(dtm), ymax(dsm))
  new_extent <- c(new_xmin, new_xmax, new_ymin, new_ymax)
  dtm <- crop(dtm, new_extent)
  dsm <- crop(dsm, new_extent)
  outputs <- list(
    dsm = dsm,         # canopy height
    dtm = dtm
  )
  res <- stack(outputs)
  return(res) 
  }


