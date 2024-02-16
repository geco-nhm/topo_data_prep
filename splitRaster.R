# split function https://stackoverflow.com/questions/29784829/r-raster-package-split-image-into-multiples
SplitRas <- function(raster,ppside,save,plot){
  h        <- ceiling(ncol(raster)/ppside)
  v        <- ceiling(nrow(raster)/ppside)
  agg      <- aggregate(raster,fact=c(h,v))
  agg[]    <- 1:ncell(agg)
  agg_poly <- rasterToPolygons(agg)
  names(agg_poly) <- "polis"
  r_list <- list()
  for(i in 1:ncell(agg)){
    e1          <- extent(agg_poly[agg_poly$polis==i,])
    r_list[[i]] <- crop(raster,e1)
  }
  if(save==T){
    for(i in 1:length(r_list)){
      writeRaster(r_list[[i]],filename=paste("SplitRas",i,sep=""),
                  format="GTiff",datatype="FLT4S",overwrite=TRUE)  
    }
  }
  if(plot==T){
    par(mfrow=c(ppside,ppside))
    for(i in 1:length(r_list)){
      plot(r_list[[i]],axes=F,legend=F,bty="n",box=FALSE)  
    }
  }
  return(r_list)
}


# crop to small extent
plot(e1_ground_mask)
testing <- crop(e1_ground_mask, oppkuven_mask_e1)
plot(testing)
testing1 <- SplitRas(raster=testing,ppside=3,save=TRUE,plot=TRUE)
plot(e1_canopy_mask)
testing <- crop(e1_canopy_mask, oppkuven_mask_e1)
plot(testing)
testing2 <- SplitRas(raster=testing,ppside=3,save=TRUE,plot=TRUE)
par(mfrow=c(1,2))
plot(testing1[[2]])
plot(testing2[[2]])

