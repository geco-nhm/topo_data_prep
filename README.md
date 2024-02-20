### Introduction
Documentation of the pipeline for processing geospatial data and deriving topographic variables, focusing on terrain and canopy data.



The file [data_prep_execute.R](data_prep_execute.R) is calling sub-functions in the pipeline process:

1. **Library Loading**: The script begins by loading necessary libraries: `raster`, `RSAGA`, `sf`, and `terra`. These are essential for handling and analyzing geospatial data.

2. **Processing Localities**: The script defines a vector `locality` with elements like "e1", "e2", "w1", "w2", and iterates over these localities in a loop.

3. **Loading Raster Data**: For each locality, it constructs paths to two raster files - a Digital Terrain Model (DTM) and a Digital Surface Model (DSM), and loads these files.

4. **Making Rasters Stackable**: The script calls a custom function [make.stackable.R](make.stackable.R) to make these rasters stackable, aligning them for further analysis.

5. **Masking with Shapefile**: The script reads a shapefile for each locality, adjusts it by dropping Z and M dimensions, and uses it to crop and mask the raster stack.

6. **Derived Data from Elevation Model using RSAGA** : The script sources a function [derive.vars.R](derive.vars.R) to derive various variables from the elevation model (DTM). This could include calculations like slope, aspect, curvature, etc.

7. **Tree Density Index** : It calculates tree density indices at different resolutions by sourcing a custom function [tree.density.R](tree.density.R). This involves running the function for multiple pixel sizes to achieve resolutions of 1m, 3m, and 9m.

8. **Resampling**: The script resamples the derived rasters for three different resolutions (1, 5, and 9 meters) and saves the resampled rasters in a specified directory [resample.raster.R](resample.raster.R).



### Functions Used in [derive.vars.R](derive.vars.R)
1. `rsaga.slope.asp.curv`
2. `rsaga.wetness.index`
3. `rsaga.hillshade`
4. `rsaga.pisr2`
5. `rsaga.topdown.processing`

### Code Structure
The code is structured to perform multiple geospatial analyses using RSAGA functions. Each segment of the code performs a specific task and saves the output in a derived data folder.

#### 1. Slope, Aspect, and Curvature Analysis
- **Function:** `rsaga.slope.asp.curv`
- **Purpose:** Derives topographical variables such as slope, aspect, and curvature from the digital elevation model (DEM).
- **Outputs:** 
  - Slope (`_slope.tif`)
    - Represents the steepness or degree of incline of the terrain. The slope is measured in degrees from the horizontal. Higher values indicate steeper terrain. Slope is a key factor in hydrology, soil erosion, vegetation distribution, and landform development.
  - Aspect (`_aspect.tif`)
    - Denotes the compass direction that a terrain surface faces. Aspect is measured in degrees from North (0°), with East at 90°, South at 180°, and West at 270°. Aspect can influence microclimate conditions like sunlight exposure and moisture retention.
  - Convergence Index (`_cgene.tif`)
    - Indicates areas of convergent and divergent flow on a landscape. Positive values typically represent convergent features like valleys or depressions, where water might accumulate. Negative values suggest divergent features like ridges or hilltops.
  - Profile Curvature (`_cprof.tif`)
    - Relates to the curvature of the terrain in the direction of the slope. It affects the acceleration and deceleration of flow across the surface, influencing erosion and deposition processes.
  - Plan Curvature (`_cplan.tif`)
    - Represents the curvature of the terrain perpendicular to the slope direction. It influences the convergence and divergence of flow across a surface, affecting soil moisture and the distribution of vegetation.
- **Units:** Degrees
- **Reference:** [RSAGA Documentation - rsaga.slope.asp.curv](https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.slope.asp.curv)

#### 2. Wetness Index Calculation
- **Function:** `rsaga.wetness.index`
- **Purpose:** Calculates the SAGA Wetness Index, which is a modification of the Topographic Wetness Index (TWI) with a more realistic approach to potential soil wetness.
- **Output:** Soil Wetness Index (`_swi.tif`)
  - An estimate of potential soil moisture content based on topography. Higher values indicate areas more likely to be wetter, such as valley bottoms. It is useful in hydrological modeling, soil science, and ecological studies.
- **Reference:** [RSAGA Documentation - rsaga.wetness.index](https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.wetness.index)

#### 3. Hillshade Generation
- **Function:** `rsaga.hillshade`
- **Purpose:** Creates a hillshade raster from the DEM using the Analytical Hillshading method.
- **Output:** Hillshade (`_hillshade.tif`)
    - A shaded relief map providing a 3D visual effect of terrain. It's created by simulating illumination over the terrain surface. Hillshade helps in visualizing and analyzing terrain features like slopes, ridges, valleys, and patterns of landform.
- **Parameters:** Standard method, Azimuth 315°, Declination 45°, Exaggeration factor 4
- **Reference:** [RSAGA Documentation - rsaga.hillshade](https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.hillshade)

#### 4. Solar Radiation Analysis (Summer Period)
- **Function:** `rsaga.pisr2`
- **Purpose:** Calculates direct, diffuse, and total insolation for the summer period.
- **Period:** June 21, 2023, to September 21, 2023
- **Period:** 6hourly intervals per day
- **Outputs:** 
  - Direct Solar Radiation (`_insolation.direct.tif`)
      - Represents the amount of solar radiation received directly from the sun, not accounting for atmospheric scattering. It is crucial for studies related to solar energy potential, ecological habitats, and climatic analysis.
  - Diffuse Solar Radiation (`_insolation.diffuse.tif`)
      - This is the solar radiation received from the sky dome without direct sunlight, scattered by the atmosphere. It's important in understanding the overall solar radiation balance, especially in shaded or cloudy areas.
  - Total Solar Radiation (`_insolation.total.tif`)
      - The sum of direct and diffuse solar radiation. It represents the total potential solar energy received over a surface, useful for solar energy projects, agricultural planning, and climate studies.
- **Reference:** [RSAGA Documentation - rsaga.pisr2](https://rdrr.io/github/r-spatial/RSAGA/man/rsaga.pisr2.html)
- **Alternative Parameters:** [GIS StackExchange Discussion](https://gis.stackexchange.com/questions/307426/rsaga-error-using-rsaga-pisr2)

#### 5. Catchment Area Calculation
- **Function:** `rsaga.topdown.processing`
- **Purpose:** Computes the local catchment area using top-down processing algorithms.
- **Output:** Catchment Area (`_catch_area.tif`)
  - Indicates the size of the area contributing runoff to each point on the terrain. Larger values denote areas that collect more water, crucial for watershed management, flood prediction, and environmental conservation.
- **Method:** Multiple Flow Direction (MFD)
- **Reference:** [RSAGA Documentation - rsaga.topdown.processing](https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.topdown.processing)


### Function [tree.density.R](tree.density.R)
this part of the script calculates canopy height, separates vegetation above 1m, calculates mean canopy height and tree density, making use of raster data operations and moving window analysis.
#### 1. Sub-functions
- **meanFunction**
  - **Description:** Calculates the mean of input values, ignoring `NA` values. This function is likely used to compute average values within a specified area or window.
- **canopyRatio**
  - **Description:** Computes the ratio of canopy cells (cells with vegetation) to total cells within a matrix. This function is integral in determining the proportion of an area covered by canopy, which is a key indicator of tree density.

#### 2. Main Function: tree_density
- **Purpose:** Calculates various metrics related to tree density based on input DEM (Digital Elevation Model), DSM (Digital Surface Model), and other parameters.
- **Process:**
  - **CHM Calculation:** Creates a Canopy Height Model by subtracting DEM from DSM, identifying vegetation height above ground.
  - **Vegetation Thresholding:** Identifies vegetation higher than 1 meter (considered as canopy).
  - **Density Calculation:**
    - Uses a moving window approach with a user-defined size.
    - Computes the mean number of canopy pixels (above 1m) within each window.
    - Calculates canopy density using the `canopyRatio` function, indicating the proportion of canopy cover within each window.
- **Outputs:** The function produces several outputs:
  - `chm`: Canopy Height Model indicating the height of vegetation above the ground.
  - `veg_1`: Binary raster indicating presence of vegetation above 1m.
  - `mean_chm`: Mean count of canopy pixels (above 1m) within the moving window across the study area.
  - `chm_density_1`: Canopy density calculated as the ratio of canopy cells to total cells in each moving window.


### File [extract.values.R](extract.values.R)
at the end, we exctract values to points at different resolutions. 
1. **Read Geospatial Data**:
   - `ina_pts` and `mari_points` are read as spatial dataframes from GeoPackage files using `st_read`. 

2. **Raster Data Processing for INA Sites**:
   - The script first handles raster data for the INA sites.
   - It reads and stacks raster layers from a specific directory (`02_derived/`) into `w2_layers`.
   - Extracts values from these raster layers at the points specified in `ina_pts`.
   - Combines (`cbind`) these values with `ina_pts`.

3. **Raster Data Processing at Different Resolutions (1m, 5m, 9m)**:
   - Iterates over different spatial resolutions (1m, 5m, 9m).
   - For each resolution, it reads and stacks corresponding raster layers from `03_resampled/` directory.
   - Extracts values from these raster layers at INA points.
   - Combines these values with the existing `ina_pts_t` dataframe.

4. **Raster Data Processing for MARI Sites**:
   - Similar steps are repeated for MARI sites with different raster layer patterns (`w1_`, `w2_`, `e1_`, `e2_`).
   - For each pattern and site combination, it reads and stacks raster layers, extracts values at MARI points, and combines these values with the MARI points data.

5. **Renaming Columns**:
   - A custom function `rename_cols` is defined to rename the column names of dataframes.
   - It removes the first `n` characters from each column name unless they are in `omit_columns`.
   - If a column name starts with any of the patterns in `patterns`, it removes the 4th, 5th, and 6th characters.

6. **Applying Renaming Function and Combining Data**:
   - The renaming function is applied to each MARI site dataframe.
   - These dataframes are then combined into a single dataframe using `rbind`.

7. **Writing Output to CSV**:
   - Finally, the combined dataframes for INA and MARI sites are written to CSV files in the `05_output` directory.

This script demonstrates advanced data manipulation involving spatial data and raster layers, tailored for ecological or geographical data processing. It efficiently handles multiple datasets, applies conditional logic for data transformation, and aggregates results into cohesive outputs.

##### this documentation was created with the help of Chat GPT
