#### CREATE DERIVED DATA FROM ELEVATION MODEL ####
# rum SAGA gis funcitons to create derived variables

#TODO make this for path creation in windows and mac
#TODO now this requires forward slash on output.dir

derive.vars <- function(output.dir, dtm, layer.code){
      env <- rsaga.env(workspace = "D:/Drone_DATA/2023_Oppkuven/",
        path = "C:/Program Files/saga-7.4.0_x64")   
      # derive topo variables in RSAGA
      rsaga.slope.asp.curv(in.dem = dtm,
                           out.slope  = paste0(output.dir, layer.code, "_slope.tif"),
                           out.aspect = paste0(output.dir, layer.code, "_aspect.tif"),
                           out.cgene =  paste0(output.dir, layer.code, "_cgene.tif"),
                           out.cprof =  paste0(output.dir, layer.code, "_cprof.tif"),
                           out.cplan =  paste0(output.dir, layer.code, "_cplan.tif"),
                           unit.slope = "degrees",
                           unit.aspect = "degrees", env = env)
      #https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.slope.asp.curv
      
      
      rsaga.wetness.index(in.dem = dtm,
                          out.wetness.index = paste0(output.dir, layer.code, "_swi.tif"),
                          env = env)
      # Details The SAGA Wetness Index is similar to the Topographic Wetness Index (TWI), but it is based on a modified catchment area calculation (out.mod.carea), which does not treat the flow as a thin film as done in the calculation of catchment areas in conventional algorithms. As a result, the SWI tends to assign a more realistic, higher potential soil wetness than the TWI to grid cells situated in valley floors with a small vertical distance to a channel.
      # https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.wetness.index
      
      a.hillshade <- rsaga.hillshade(in.dem = dtm,
                                     out.grid = paste0(output.dir, layer.code, "_hillshade.tif"),
                                     method = "standard",
                                     azimuth = 315,
                                     declination = 45,
                                     exaggeration = 4,
                                     env = env)
      #The Analytical Hillshading algorithm is based on the angle between the surface and the incoming light beams, measured in radians.
      #https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.hillshade
      
      # rsaga.insolation(in.dtm = area.name,
      #                  out.direct = paste0("SAGA/", layer.code, "_insolation.direct.tif"),
      #                  out.diffuse = paste0("SAGA/", layer.code, "_insolation.diffuse.tif"),
      #                  out.total = paste0("SAGA/", layer.code, "_insolation.total.tif"),
      #                  time.step = 4, day.step = 1, days = 172:264, 
      #                  env = env)
      
      #in the Northern Hemisphere, summer might be considered from around June 21st (day 172) to September 21st (day 264).
      #https://www.rdocumentation.org/packages/RSAGA/versions/1.4.0/topics/rsaga.insolation
      # Calculation of incoming solar radiation (insolation). Based on the SADO (System for the Analysis of Discrete Surfaces) routines developed by Boehner & Trachinow.
      
      rsaga.pisr2(in.dem = dtm,
                  out.direct.grid = paste0(output.dir, layer.code, "_insolation.direct.tif"),
                  out.diffuse.grid = paste0(output.dir, layer.code, "_insolation.diffuse.tif"),
                  out.total.grid = paste0(output.dir, layer.code, "_insolation.total.tif"),
                  out.duration =  paste0(output.dir, layer.code, "_insolation.duration.tif"),
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
      rsaga.topdown.processing(in.dem = dtm,
                               out.carea = paste0(output.dir, layer.code, "_catch_area.tif"),
                               method = "mfd", 
                               env = env)
}