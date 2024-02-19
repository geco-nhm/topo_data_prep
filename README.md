To create documentation for the provided code, you'll want to organize and present the information in a way that is clear, informative, and easy for other users or future you to understand. Here's a structured approach to documenting this code:

### Title
**Documentation for Derived Data Creation from Elevation Model Using RSAGA**

### Introduction
This document provides detailed information about the process of creating derived data from an elevation model using various RSAGA functions. The RSAGA package is an interface between R and the SAGA GIS software.

### Functions Used
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
  - Aspect (`_aspect.tif`)
  - Convergence Index (`_cgene.tif`)
  - Profile Curvature (`_cprof.tif`)
  - Plan Curvature (`_cplan.tif`)
- **Units:** Degrees
- **Reference:** [RSAGA Documentation - rsaga.slope.asp.curv](https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.slope.asp.curv)

#### 2. Wetness Index Calculation
- **Function:** `rsaga.wetness.index`
- **Purpose:** Calculates the SAGA Wetness Index, which is a modification of the Topographic Wetness Index (TWI) with a more realistic approach to potential soil wetness.
- **Output:** Soil Wetness Index (`_swi.tif`)
- **Reference:** [RSAGA Documentation - rsaga.wetness.index](https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.wetness.index)

#### 3. Hillshade Generation
- **Function:** `rsaga.hillshade`
- **Purpose:** Creates a hillshade raster from the DEM using the Analytical Hillshading method.
- **Output:** Hillshade (`_hillshade.tif`)
- **Parameters:** Standard method, Azimuth 315°, Declination 45°, Exaggeration factor 4
- **Reference:** [RSAGA Documentation - rsaga.hillshade](https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.hillshade)

#### 4. Solar Radiation Analysis (Summer Period)
- **Function:** `rsaga.pisr2`
- **Purpose:** Calculates direct, diffuse, and total insolation for the summer period.
- **Period:** June 21, 2023, to September 21, 2023
- **Period:** 6hourly intervals per day
- **Outputs:** 
  - Direct Solar Radiation (`_insolation.direct.tif`)
  - Diffuse Solar Radiation (`_insolation.diffuse.tif`)
  - Total Solar Radiation (`_insolation.total.tif`)
  - Duration of Insolation (`_insolation.duration.tif`)
- **Reference:** [RSAGA Documentation - rsaga.pisr2](https://rdrr.io/github/r-spatial/RSAGA/man/rsaga.pisr2.html)
- **Alternative Parameters:** [GIS StackExchange Discussion](https://gis.stackexchange.com/questions/307426/rsaga-error-using-rsaga-pisr2)

#### 5. Catchment Area Calculation
- **Function:** `rsaga.topdown.processing`
- **Purpose:** Computes the local catchment area using top-down processing algorithms.
- **Output:** Catchment Area (`_catch_area.tif`)
- **Method:** Multiple Flow Direction (MFD)
- **Reference:** [RSAGA Documentation - rsaga.topdown.processing](https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.topdown.processing)

### General Notes
- Ensure all necessary packages (`RSAGA`, `raster`, etc.) are installed and loaded in R.
- Replace `file.name` and `layer.code` with the specific file names and layer codes relevant to your project.
- The output files are saved in a folder named `03_derived` with a naming convention based on the `layer.code`.

### Conclusion
This documentation provides an overview of the script used for deriving various geospatial data from an elevation model using RSAGA. It is designed to guide the user through each step of the process, explaining the purpose, function usage, and output of each part of the script.



Outputs Description
1. Slope (_slope.tif)
Description: Represents the steepness or degree of incline of the terrain. The slope is measured in degrees from the horizontal. Higher values indicate steeper terrain. Slope is a key factor in hydrology, soil erosion, vegetation distribution, and landform development.
2. Aspect (_aspect.tif)
Description: Denotes the compass direction that a terrain surface faces. Aspect is measured in degrees from North (0°), with East at 90°, South at 180°, and West at 270°. Aspect can influence microclimate conditions like sunlight exposure and moisture retention.
3. Convergence Index (_cgene.tif)
Description: Indicates areas of convergent and divergent flow on a landscape. Positive values typically represent convergent features like valleys or depressions, where water might accumulate. Negative values suggest divergent features like ridges or hilltops.
4. Profile Curvature (_cprof.tif)
Description: Relates to the curvature of the terrain in the direction of the slope. It affects the acceleration and deceleration of flow across the surface, influencing erosion and deposition processes.
5. Plan Curvature (_cplan.tif)
Description: Represents the curvature of the terrain perpendicular to the slope direction. It influences the convergence and divergence of flow across a surface, affecting soil moisture and the distribution of vegetation.
6. Soil Wetness Index (_swi.tif)
Description: An estimate of potential soil moisture content based on topography. Higher values indicate areas more likely to be wetter, such as valley bottoms. It is useful in hydrological modeling, soil science, and ecological studies.
7. Hillshade (_hillshade.tif)
Description: A shaded relief map providing a 3D visual effect of terrain. It's created by simulating illumination over the terrain surface. Hillshade helps in visualizing and analyzing terrain features like slopes, ridges, valleys, and patterns of landform.
8. Direct Solar Radiation (_insolation.direct.tif)
Description: Represents the amount of solar radiation received directly from the sun, not accounting for atmospheric scattering. It is crucial for studies related to solar energy potential, ecological habitats, and climatic analysis.
9. Diffuse Solar Radiation (_insolation.diffuse.tif)
Description: This is the solar radiation received from the sky dome without direct sunlight, scattered by the atmosphere. It's important in understanding the overall solar radiation balance, especially in shaded or cloudy areas.
10. Total Solar Radiation (_insolation.total.tif)
Description: The sum of direct and diffuse solar radiation. It represents the total potential solar energy received over a surface, useful for solar energy projects, agricultural planning, and climate studies.
11. Duration of Insolation (_insolation.duration.tif)
Description: Shows the duration for which each part of the terrain is exposed to solar radiation. It can be used for analyzing variations in solar exposure due to terrain features.
12. Catchment Area (_catch_area.tif)
Description: Indicates the size of the area contributing runoff to each point on the terrain. Larger values denote areas that collect more water, crucial for watershed management, flood prediction, and environmental conservation.


### File extract.values.R 
1. Read Geospatial Data:
* ina_pts and mari_points are read as spatial dataframes from GeoPackage files using st_read.
* Raster Data Processing for INA Sites:
* The script first handles raster data for the INA sites.
* It reads and stacks raster layers from a specific directory (02_derived/) into w2_layers.
* Extracts values from these raster layers at the points specified in ina_pts.
* Combines (cbind) these values with ina_pts.
2. Raster Data Processing at Different Resolutions (1m, 5m, 9m):

* Iterates over different spatial resolutions (1m, 5m, 9m).
* For each resolution, it reads and stacks corresponding raster layers from 03_resampled/ directory.
* Extracts values from these raster layers at INA points.
* Combines these values with the existing ina_pts_t dataframe.

3. Raster Data Processing for MARI Sites:
* Similar steps are repeated for MARI sites with different raster layer patterns (w1_, w2_, e1_, e2_).
* For each pattern and site combination, it reads and stacks raster layers, extracts values at MARI points, and combines these values with the MARI points data.
4. Renaming Columns:
A custom function rename_cols is defined to rename the column names of dataframes.
It removes the first n characters from each column name unless they are in omit_columns.
If a column name starts with any of the patterns in patterns, it removes the 4th, 5th, and 6th characters.
Applying Renaming Function and Combining Data:

The renaming function is applied to each MARI site dataframe.
These dataframes are then combined into a single dataframe using rbind.
Writing Output to CSV:

Finally, the combined dataframes for INA and MARI sites are written to CSV files in the 05_output directory.