library(stringr)
library(pbapply)

#### EXTRACT POINTS ####

# at each site we extract first the original 0.1m resolution data *these are without prefix
# then we loop through resolutions of 1,5,9m and extract for each var *these are with prefix
# then we fix names so that site code is omitted 
# and last we combine into DF and save as csv and geopackage


DataDir <- "D:/Drone_DATA/2023_Oppkuven/"

mari_points <- st_read(file.path(DataDir, "Mari_sites.gpkg"))
mari_pts <- mari_points
mari_pts$site <- str_replace_all(mari_pts$site, c("2" = "w2", "1" = "w1", "3" = "e1", "4" = "e2"))






SiteResExtraction <- function(Sites = tolower(c("W1", "W2", "E1", "E2")),
                              Resolutions = c(1,5,9),
                              locations = mari_pts # An sf point object. must contain a column called "site" with entries corresponding to the Sites argument. Must also contain a plot_id column
                              ){
  Extraction_ls <- pblapply(Sites, FUN = function(Site_iter){
    ### Base Resolution ----
    layer_list <- list.files(paste0(DataDir,'/02_derived/'), pattern=paste0(Site_iter, "_"), full.names=TRUE)
    # read in raster layers and stack them into raster stack
    layers <- stack(layer_list)
    # Extract values
    mari_pts_values <- cbind(subset(locations, site==Site_iter)[,c("plot_id", "site")], extract(layers, subset(locations, site==Site_iter)))
    mari_pts_values <- data.frame(mari_pts_values)
    mari_pts_values <- mari_pts_values[,-ncol(mari_pts_values)]
    colnames(mari_pts_values) <- str_replace(colnames(mari_pts_values), pattern = paste0(Site_iter, "_"), replacement = "")
    
    ### Resampled Resolutions ----
    Res_ls <- lapply(Resolutions, FUN = function(Res_iter){
      layer_list_res <- list.files(paste0(DataDir,'/03_resampled/'), pattern=paste0(Res_iter,"m_",Site_iter, "_"), full.names=TRUE)
      # read in raster layers and stack them into raster stack
      layers_res <- stack(layer_list_res)
      mari_res_pts_values <- cbind(subset(locations, site==Site_iter)[,c("plot_id", "site")], extract(layers_res, subset(locations, site==Site_iter)))
      mari_res_pts_values <- data.frame(mari_res_pts_values)
      mari_res_pts_values <- mari_res_pts_values[,-ncol(mari_res_pts_values)]
      colnames(mari_res_pts_values) <- str_replace(colnames(mari_res_pts_values), pattern = paste0(Site_iter, "_"), replacement = "")
      mari_res_pts_values
    })
    if(length(Resolutions) > 1){ # if the user queried more than one resolution
      Res_df <- Reduce(function(x, y) merge(x, y), Res_ls) # we merge all the data frames in the resampled extraction list (they are merged by site and plot_id)
    }else{ # if the user specified just one resolution
      Res_df <- Res_ls[[1]] # we just use that resampling resolution extraction (there won't be anything to merge)
    }
    
    ### Return objects to higher level list ----
    list(Base = mari_pts_values,
         Resampled = Res_df)
  })
  Base_df <- do.call(rbind, lapply(Extraction_ls, "[[", "Base"))
  Res_df <- do.call(rbind, lapply(Extraction_ls, "[[", "Resampled"))
  Return_df <- merge(Base_df, Res_df)
  return(Return_df)
}

test <- SiteResExtraction(Sites = "w2")
