# Foundation / source of code: https://r-video-tutorial.blogspot.com/2015/08/spatio-temporal-kriging-in-r.html
library(gstat)
library(sp)
library(spacetime)
library(raster)
library(rgdal)
library(rgeos) 
library(RPostgres)
library(sf)
library(future)
library(future.apply)

# Set working directory
setwd("R:/Analysis/Scripts/rensvwijk/R_scripts/data")

# Set parameters of script
VARIABLE <- "temperature_mean"
START_DATE <- "2023-01-01"
END_DATE <- "2023-01-14"
START_TIME <- "00:00"
END_TIME <-  "23:00"
K_STATIONS <- 50
BOUNDARY_DATA = "TerraClim_2020_v04.shp"

# Create Folder
FOLDER_NAME <- gsub("[-:]", "", sprintf("%s/Model_%s_%s_%s_%s", Sys.Date (),START_DATE, START_TIME,END_DATE, END_TIME))
FOLDER <- file.path("R:/Analysis/Scripts/rensvwijk/R_scripts/Space-Time Kriging/Results", FOLDER_NAME)

if (!file.exists(FOLDER)) {
  dir.create(FOLDER, recursive = TRUE) 
}else {
  cat("WARNING: The folder", FOLDER, "already exists.\n")
  cat("Do you want to continue and overwrite the existing folder'? (y/n): ")
  
  user_input <- readline(prompt = "")
  
  if (tolower(user_input) == "y") {
    # Remove the existing folder
    unlink(FOLDER, recursive = TRUE)
    
    # Recreate the folder
    dir.create(FOLDER, recursive = TRUE)
    cat("The folder", FOLDER, "has been overwritten.\n")
    cat("Operation continues.")
  } else {
    stop("Operation canceled. No changes were made to the existing folder.\n")
  }
}

######################

# DATABASE DATA 
# Set database information 
db_host <- "xxx.xxx.xx.xxx"
db_port <- 5432
db_name <- "XXXXX"
db_user <- "XXXXX"
db_password <- "XXXXX"

# Establish a connection
con <- dbConnect(
  RPostgres::Postgres(),
  host = db_host,
  port = db_port,
  dbname = db_name,
  user = db_user,
  password = db_password 
)

# Query data from database
query <- sprintf('SELECT date, time, %s, "Weather_Stations_v3".station_id, "Weather_Stations_v3".latitude, "Weather_Stations_v3".longitude, public.ST_AsText(geom) AS geom_text FROM "Hourly_Clean_v3" INNER JOIN "Weather_Stations_v3" ON "Hourly_Clean_v3".station_id = "Weather_Stations_v3".station_id WHERE "date" >= \'%s\'::date AND "date" <= \'%s\'::date AND "covered" = \'NO\'', VARIABLE, START_DATE, END_DATE)

# Execute the query and fetch the results
dbListTables(con)
data <- dbGetQuery(con, query)

######################

# DATA PREPROCESSING
# Standardize date/time data to UNIX time 
data$time <- as.numeric(data$time)
data$TIME <- as.POSIXlt(paste(data$date, sprintf("%02d:%02d", data$time %/% 100, data$time %% 100), "Africa/Johannesburg"), origin="1970-01-01")

# Rename temperature column
data$TEMP <- data$temperature_mean

# Obtain coordinates and transform to manageable format
data$LAT <- as.numeric(data$lat)
data$LON <- as.numeric(data$long)

# Exclude NoData from dataset
data <- na.omit(data) 

# Choose two time periods to minimize size of dataset for kriging
subset <- data[data$TIME >= as.POSIXct(sprintf("%s %s SAST", START_DATE, START_TIME)) & data$TIME <= as.POSIXct(sprintf("%s %s SAST", END_DATE, END_TIME)), ]

# Create a SpatialPointsDataFrame
coordinates(subset)=~LON+LAT
projection(subset)=CRS("+init=epsg:4326")

#Transform into Mercator Projection
temp.UTM <- spTransform(subset,CRS("+init=epsg:3395")) 

# Change projection of boundary data from original to UTM 
area <- shapefile(BOUNDARY_DATA)
area.UTM <- spTransform(area,CRS("+init=epsg:3395")) 

# Crop file to research area
temp.UTM.cropped <- crop(temp.UTM, extent(area.UTM)) 

# Create SpatialPoints object with locations of sensors at any given time
tempSP <- SpatialPoints(temp.UTM.cropped@coords,CRS("+init=epsg:3395")) 

# Remove duplicate points (data)
tempDF <- data.frame(TEMP=temp.UTM.cropped$temperature_mean) 

# Remove duplicate points (time)
tempTM <- as.POSIXct(temp.UTM.cropped$TIME,tz="Africa/Johannesburg") 

# Combine original, and cleaned data (data and time)
timeDF <- STIDF(tempSP,tempTM,data=tempDF) 

######################

# PLOTTING OF DATA
# Create plot saver
plot_and_save <- function(plot, filename) {
  filename <- sprintf("%s.png", filename)
  png(filename = file.path(FOLDER, filename))
  print(plot)
  dev.off()
  print(plot)
}

# Plot interim dataset
plot_and_save(stplot(timeDF), "Stations")

# Create Shapefile of Western Cape province as SpatialPolygons class
WC_shp <- st_read("Province2006.shp") 
WC_polygon <- WC_shp[WC_shp$code == "WC", ]
WC_polygon.UTM <- st_transform(WC_polygon, crs = "+init=epsg:3395")
WC_spatialpolygons <- as(WC_polygon.UTM, "Spatial")

# Create one big area for TerraClim regions
total_area.UTM <- gUnaryUnion(area.UTM)

# Intersection of two polygons
WC_TerraClim <- gIntersection(WC_spatialpolygons, total_area.UTM)

# Intersect SpatialPoints with TerraClim region
WC_tempSP <- gIntersection(tempSP, WC_TerraClim)

# Plot map of South Africa
plot_TerraClim <- function() {
  # Create Shapefile of South Africa as SpatialPolygon class
  SA_shp.UTM <- st_transform(WC_shp, crs = "+init=epsg:3395")
  SA_spatialpolygons <- as(SA_shp.UTM, "Spatial")
  
  # Plot South Africa on map with WC Province and TerraClim Regions
  plot(SA_spatialpolygons, main = "Locations of the TerraClim region in South Africa", col="lightgrey")
  plot(WC_spatialpolygons, col = "darkgrey", add=T)
  plot(WC_TerraClim, col = "lightgreen", add=T)
  
  # Add black outline and legend
  box()
  legend(
    "topleft",
    legend = c("TerraClim Region", "Western Cape Province", "South Africa"),
    col = c("lightgreen", "darkgrey", "lightgrey"),
    pch = c(15, 15, 15),  # Use pch 15 for squares in the legend
    title = "Legend"
  )
}

# Plot and Save 
plot_and_save(plot_TerraClim(),"TerraClimRegion_WC")

# Plot weather stations on final Map
plot_weather_stations <- function() {
  # Plot Weather Stations on map of WC Province and TerraClim Regions
  plot(WC_spatialpolygons, main = "Locations of weather stations in the TerraClim region", col="grey")
  plot(WC_TerraClim, col = "lightgreen", add=T)
  plot(WC_tempSP, pch = 20, col = "black", add = T)
  
  # Add black outline and legend
  box()
  legend(
    "topright",
    legend = c("Weather Station", "TerraClim Region", "Western Cape Province"),
    col = c("black", "lightgreen", "lightgrey"),
    pch = c(20, 15, 15),  # Use pch 15 for squares in the legend
    title = "Legend"
  )
}

# Plot and Save 
plot_and_save(plot_weather_stations(),"WeatherStations_WC")

######################

# VARIOGRAM
# Compute variogram
var <- variogramST(
  TEMP~1,
  data=timeDF,
  tunit="hours",
  twindow=3,
  cores=24
)

# Plot variograms
plot_and_save(plot(var, wireframe=T, scales=list(arrows=F), zlab = "temperature"), "Variogram")

# Plot variograms (normal type)
plot_and_save(plot(var,wireframe=T, zlab = "temperature"), "Variogram_normal") 

# Plot variograms (map type)
plot_and_save(plot(var,map=T, zlab = "temperature"), "Variogram_map") 

# Plot variograms (3D wireframe type)
plot_and_save(plot(var,wireframe=T, zlab = "temperature"), "Variogram_3D") 

######################

# VARIOGRAM MODELLING
# Five options for variogram modelling:
## 1 - Separable 
## 2 - Product sum
## 3 - Metric
## 4 - Sum metric
## 5 - Simple sum metric

# OPTION 1 - SEPARABLE
# Provide spatial component, temporal component and sill
separable <- vgmST(
  "separable", 
  space = vgm(-60,"Sph", 500, 1),
  time = vgm(35,"Sph", 500, 1), 
  sill=0.56
) 

plot_and_save(plot(var,separable,map=F,y="temperature"), "Vgm_separable")

# Model Fit
separable_Vgm <- fit.StVariogram(
  var, 
  separable, 
  fit.method=0
)
sprintf("Option 1 - separable model fit = %f", attr(separable_Vgm,"MSE"))

# Check Model Fit (alternative option)
separable_alt_Vgm <- fit.StVariogram(
  var, 
  separable, 
  fit.method=11,
  method="L-BFGS-B", 
  stAni=5, 
  )
sprintf("Option 1 (alternative) - separable model fit = %f", attr(separable_alt_Vgm,"MSE"))

# Check parameters of model
print(extractPar(separable_Vgm))

# OPTION 2 - PRODUCT SUM
# Provide spatial component, temporal component and value of parameter k
prodSumModel <- vgmST(
  "productSum",
  space = vgm(1, "Exp", 150, 0.5),
  time = vgm(1, "Exp", 5, 0.5),
  k = 50
) 

# Model Fit
prodSumModel_Vgm <- fit.StVariogram(
  var, 
  prodSumModel,
  method = "L-BFGS-B"
  )
sprintf("Option 2 - product sum model fit = %f", attr(prodSumModel_Vgm,"MSE"))

# Plot variogram
plot_and_save(plot(var,prodSumModel_Vgm,map=F,y="temperature"), "Vgm_prodSumModel") 

# Check parameters of model
print(extractPar(prodSumModel_Vgm))

# OPTION 3 - METRIC
# Provide value of parameter k
metric <- vgmST(
  "metric", 
  joint = vgm(50,"Mat", 500, 0), 
  stAni=200
) 

# Model Fit
metric_Vgm <- fit.StVariogram(
  var, 
  metric, 
  method="L-BFGS-B"
)
sprintf("Option 3 - metric model fit = %f", attr(metric_Vgm,"MSE"))

# Plot variogram
plot_and_save(plot(var,metric_Vgm,map=F,y="temperature"), "Vgm_metric") 

# Check parameters of model
print(extractPar(metric_Vgm))

# OPTION 4 - SUM METRIC
# Provide all parameters independently
sumMetric <- vgmST(
  "sumMetric", 
  space = vgm(psill=5,"Sph", range=500, nugget=0),
  time = vgm(psill=500,"Sph", range=500, nugget=0), 
  joint = vgm(psill=1,"Sph", range=500, nugget=10), 
  stAni=500
) 

# Model Fit
sumMetric_Vgm <- fit.StVariogram(
  var, 
  sumMetric, 
  method="L-BFGS-B",
  tunit="hours"
)
sprintf("Option 4 - sum metric model fit = %f", attr(sumMetric_Vgm,"MSE"))

# Plot variogram
plot_and_save(plot(var,sumMetric_Vgm,map=F,y="temperature"), "Vgm_sumMetric") 

# Check parameters of model
print(extractPar(sumMetric_Vgm))

# OPTION 5 - SIMPLE SUM METRIC
# Provide all parameters independently
SimplesumMetric <- vgmST(
  "simpleSumMetric",
  space = vgm(5,"Sph", 500, 0),
  time = vgm(500,"Sph", 500, 0), 
  joint = vgm(1,"Sph", 500, 0), 
  nugget=1, 
  stAni=500
) 

# Model Fit
SimplesumMetric_Vgm <- fit.StVariogram(
  var, 
  SimplesumMetric,method = "L-BFGS-B"
)
sprintf("Option 5 - simple sum metric model fit = %f", attr(SimplesumMetric_Vgm,"MSE"))

# Plot variogram
plot_and_save(plot(var,SimplesumMetric_Vgm,map=F,y="temperature"), "Vgm_simpleSumMetric") 

# Check parameters of model
print(extractPar(SimplesumMetric_Vgm))

######################

# COMPARE MODELS
plot_and_save(plot(var,list(separable_Vgm, prodSumModel_Vgm, metric_Vgm, sumMetric_Vgm, SimplesumMetric_Vgm),all=T,wireframe=T,xlab="", ylab="",zlab=""), "Variogram_Compare")  

# Calculate MSE per model and store data frame as CSV
model_names <- c("separable_Vgm", "prodSumModel_Vgm", "metric_Vgm", "sumMetric_Vgm", "SimplesumMetric_Vgm")
result_df <- data.frame(
  VariableName = model_names,
  MSE = unlist(lapply(model_names, function(model) attr(get(model), "MSE")))
)
write.csv(result_df, file = sprintf("%s/modelMSE.csv", FOLDER), row.names = FALSE)
print(result_df)

# Chose best performing model
best_model <- result_df$VariableName[which.min(result_df$MSE)]
print(best_model)

######################

# COVARIATES 
# Function to add covariate from raster [FUTURE PURPOSES: THIS FUNCTION ENABLES THE ADDING OF COVARIATES IN timeDF]
addCovariateFromRaster <- function(rasterFile, covariateName) {
  # Load and process the raster
  covariateRaster <- raster(rasterFile)
  covariateRasterUTM <- projectRaster(covariateRaster, crs = CRS("+init=epsg:3395"))
  
  # Extract covariate values at station locations
  covariateValues <- extract(covariateRasterUTM, tempSP)
  
  # Add the covariate to the data frame
  timeDF[[covariateName]] <- covariateValues
  
  return(timeDF)
}

# Apply function for covariates [FUTURE PURPOSES: THESE FUNCTIONS CALLS THE FUNCTION TO ADD COVARIATES IN timeDF]
timeDF <- addCovariateFromRaster("SUDEM_40m_WC_WGS.tif", "SUDEM")
timeDF <- addCovariateFromRaster("DtC_40m_WC_WGS.tif", "DtC")
timeDF <- addCovariateFromRaster("SR_40m_WC_WGS.tif", "SR")

######################

# SPATIAL-TEMPORAL PREDICTION GRID 
# Create random grid of points around temp points
sp.grid.UTM <- spsample(
  WC_TerraClim, 
  n=10000,
  type="random"
) 

# Create function to calculate total hours [FUTURE PURPOSES: THIS COULD BE USED TO CALCULATE INTERPOLATION FOR EACH HOUR IN GIVEN TIME FRAME]
totalTimeDifference <- function() { 
  # Combine date and time strings into POSIXct date-time objects
  start_datetime <- as.POSIXct(paste(START_DATE, START_TIME), tz = "UTC")
  end_datetime <- as.POSIXct(paste(END_DATE, END_TIME), tz = "UTC")
  
  # Calculate differences in days and hours
  diff_hours <- as.numeric(difftime(end_datetime, start_datetime, units = "hours"))
  total_difference <- diff_hours + 1
  
  return(total_difference)
}

# Create a vector of Date/Times for adding temporal component for each hours [FUTURE PURPOSES: THIS COULD BE USED TO CALCULATE INTERPOLATION FOR EACH HOUR IN GIVEN TIME FRAME]
#tm.grid <- seq(as.POSIXct(sprintf("%s %s SAST", START_DATE, START_TIME)),as.POSIXct(sprintf("%s %s SAST", END_DATE, END_TIME)),length.out=totalTimeDifference()) 

# Create a vector of Date/Times for adding temporal component
tm.grid <- seq(as.POSIXct(sprintf("%s %s SAST", START_DATE, START_TIME)),as.POSIXct(sprintf("%s %s SAST", END_DATE, END_TIME)),length.out=6) 

# Merge spatial and temporal data into a spatio-temporal grid
grid.ST <- STF(sp.grid.UTM,tm.grid) 

######################

# SPACE-TIME KRIGING INTERPOLATION
# Perform space-time kriging interpolation
pred <- krigeST(
  TEMP~1, 
  #TEMP~1 + SUDEM + DtC + SR, # [FUTURE PURPOSES: THIS LINE COULD BE USED TO ADD COVARIATES INSTEAD OF ONE ABOVE]
  data=timeDF,
  modelList=get(best_model), 
  newdata=grid.ST, 
  nmax=K_STATIONS,
  main=sprintf("Prediction for: %s %s %s %s\nModel: %s", START_DATE, START_TIME,END_DATE, END_TIME, best_model)
) 

# Plot and save results
plot_and_save(stplot(pred), "krigeST") 
