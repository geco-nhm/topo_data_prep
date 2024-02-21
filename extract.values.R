#### EXTRACT POINTS ####

# at each site we extract first the original 0.1m resolution data *these are without prefix
# then we loop through resolutions of 1,5,9m and extract for each var *these are with prefix
# then we fix names so that site code is omitted 
# and last we combine into DF and save as csv and geopackage



#### INA ####
ina_pts <- st_read("Ina_sites.gpkg")
layer.code <- "w2"
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
st_write(ina_pts_t, paste0(getwd(),'/05_output/ina_pts.gpkg'),
         layer = "mari_derived_vars")


#### MARI ####
mari_points <- st_read("Mari_sites.gpkg")
mari_pts <- mari_points


#w1 ####
w1_layer_list <- list.files(paste0(getwd(),'/02_derived/'), pattern="w1_", full.names=TRUE)
# read in raster layers and stack them into raster stack
w1_layers <- stack(w1_layer_list)

# Extract values
mari_pts_values_w1 <- extract(w1_layers, subset(mari_pts, site=="1"))
# bind columns to the correct site
mari_pts_1 <- cbind(subset(mari_pts, site=="1"), mari_pts_values_w1)

#then take the three resolutions at 1,5,9m
res <- c("1m","5m","9m")
for(x in res){
  w1_layer_list_res <- list.files(paste0(getwd(),'/03_resampled/'), pattern=paste0(x,"_","w1", "_"), full.names=TRUE)
  # read in raster layers and stack them into raster stack
  w1_layers_res <- stack(w1_layer_list_res)
  
  mari_w1_res_pts_values <- extract(w1_layers_res, subset(mari_pts, site=="1"))
  mari_pts_1 <- cbind(subset(mari_pts_1, site=="1"), mari_w1_res_pts_values)
}



#w2 ####
w2_layer_list <- list.files(paste0(getwd(),'/02_derived/'), pattern="w2_", full.names=TRUE)
# read in raster layers and stack them into raster stack
w2_layers <- stack(w2_layer_list)

# Extract values
mari_pts_values_w2 <- extract(w2_layers, subset(mari_pts, site=="2"))
# bind columns to the correct site
mari_pts_2 <- cbind(subset(mari_pts, site=="2"), mari_pts_values_w2)

#then take the three resolutions at 1,5,9m
res <- c("1m","5m","9m")
for(x in res){
  w2_layer_list_res <- list.files(paste0(getwd(),'/03_resampled/'), pattern=paste0(x,"_","w2", "_"), full.names=TRUE)
  # read in raster layers and stack them into raster stack
  w2_layers_res <- stack(w2_layer_list_res)
  
  mari_w2_res_pts_values <- extract(w2_layers_res, subset(mari_pts, site=="2"))
  mari_pts_2 <- cbind(subset(mari_pts_2, site=="2"), mari_w2_res_pts_values)
}

### e1 ####
e1_layer_list <- list.files(paste0(getwd(),'/02_derived/'), pattern="e1_", full.names=TRUE)
# read in raster layers and stack them into raster stack
e1_layers <- stack(e1_layer_list)
# Extract values
mari_pts_values_e1 <- extract(e1_layers, subset(mari_pts, site=="3"))
# bind columns to the correct site
mari_pts_3 <- cbind(subset(mari_pts, site=="3"), mari_pts_values_e1)

#then take the three resolutions at 1,5,9m
res <- c("1m","5m","9m")
for(x in res){
  e1_layer_list_res <- list.files(paste0(getwd(),'/03_resampled/'), pattern=paste0(x,"_","e1", "_"), full.names=TRUE)
  # read in raster layers and stack them into raster stack
  e1_layers_res <- stack(e1_layer_list_res)
  
  mari_e1_res_pts_values <- extract(e1_layers_res, subset(mari_pts, site=="3"))
  mari_pts_3 <- cbind(subset(mari_pts_3, site=="3"), mari_e1_res_pts_values)
}



####e2 ####
e2_layer_list <- list.files(paste0(getwd(),'/02_derived/'), pattern="e2_", full.names=TRUE)
# read in raster layers and stack them into raster stack
e2_layers <- stack(e2_layer_list)

# Extract values
mari_pts_values_e2 <- extract(e2_layers, subset(mari_pts, site=="4"))
# bind columns to the correct site
mari_pts_4 <- cbind(subset(mari_pts, site=="4"), mari_pts_values_e2)

#then take the three resolutions at 1,5,9m
res <- c("1m","5m","9m")
for(x in res){
  e2_layer_list_res <- list.files(paste0(getwd(),'/03_resampled/'), pattern=paste0(x,"_","e2", "_"), full.names=TRUE)
  # read in raster layers and stack them into raster stack
  e2_layers_res <- stack(e2_layer_list_res)
  
  mari_e2_res_pts_values <- extract(e2_layers_res, subset(mari_pts, site=="4"))
  mari_pts_4 <- cbind(subset(mari_pts_4, site=="4"), mari_e2_res_pts_values)
}



# rename and combine ####
# remove the first n characters from each column name
n <- 3  # Number of characters to remove
omit_columns <- c("geometry", "site", "plot_id")
patterns <- c("X1m_", "X5m_", "X9m_") # for grepl search 
#In this example, n + 1 is the starting position from where you want to keep the characters (since you're removing the first n characters), and nchar(colnames(df)) gives the length of each column name, ensuring you keep everything from the n+1th character to the end.
# Loop through column names but omit
mari.sites <- list(mari_pts_1,mari_pts_2,mari_pts_3,mari_pts_4)

# Function to rename columns with additional condition
rename_cols <- function(df) {
    colnames(df) <- sapply(colnames(df), function(col_name) {
      if (!col_name %in% omit_columns) {
        # Check if column name starts with "X1"
        if (grepl(paste(patterns, collapse="|"), col_name)) {
          # Remove 4th, 5th, and 6th characters
          part1 <- substring(col_name, 1, 3)
          part2 <- substring(col_name, 7, nchar(col_name))
          return(paste0(part1, part2))
        } else {
          # Otherwise, remove the first n characters
          return(substring(col_name, n + 1))
        }
      } else {
        # Keep the original name for omitted columns
        return(col_name)
      }
    })
    return(df)
  }

# Apply the function to each dataframe in the list
mari.sites <- lapply(mari.sites, rename_cols)
# Perform rbind on all dataframes in the list
mari_pts_t <- do.call(rbind, mari.sites)
write.csv2(mari_pts_t, paste0(getwd(),'/05_output/mari_pts.csv'))
st_write(mari_pts_t, paste0(getwd(),'/05_output/mari_pts.gpkg'),
         layer = "mari_derived_vars")

#save column names and add variable descriptions
var_description <- colnames(mari_pts_t)
write.csv2(var_description, paste0(getwd(),'/05_output/var_description.csv'))
