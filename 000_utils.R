#
# Utilities needed for preprocessing & analysing NCDC weather dataset
# (Can be regularly updated)
#

# String for latitude-longitude projection (can be passed into CRS)
g_latlongProj_str <- "+proj=longlat +datum=WGS84"

# Function to remove temporary objects
## Convention: name starting with 'c_'
rmMyTemp_fn <- function(pattern="^c_.*$") {
  rm(list=apropos(pattern), pos=globalenv()) 
  # CAREFUL: accesses global environment. Know what you are doing
}

# Function to parse string -> POSIX
parserDate_fn <- function(str, tz) {
  if (!("package:lubridate" %in% search())) require(lubridate)
  ## PRE: string containing information about date and time (well-formatted)
  ## tz is a string for time zone -- here, restricted to US time zones.
  ## POST: UTC converted, POSIX vector
    time <- as.POSIXct(str, tz=tz)
    as_datetime(time, tz="UTC")
}

# Image saving function (png format, fixed resolution)
# (For ggplot use ggsave instead)
savePlot_fn <- function(filename, plotfn, ...) {
  ## PRE: plotfn can be plot. method or external plot functions 
  ## ... passed onto plotfn as arguments
  resultPlot <- plotfn(...)
  png(filename, width=1080, height=720)
  print(resultPlot)
  dev.off()
}

# Function to display pairwise ggplots 
# that do not necessarily share the same variables
# (May be generalized to include multiple, distinct ggplots)
pairGG_fn <- function(g1, g2, ...) {
  ## ... passed to ggmatrix fn
  plotList <- list(g1, g2)
  resultPlotMat <- ggmatrix(plotList, nrow=1, ncol=2, ...) 
  # Parameters: label / legend / aes, etc.
  resultPlotMat
}
#---
# Function that creates a binary weights matrix for time series (see Griffith, 2013)
binAdjTime_fn <- function(n, lags=NULL) {
  # only parameter is n, the dimension (number of time points)
  adjMat <- matrix(0, n, n)
  if (is.null(lags)) {
    for (i in 1:n) {
      if (i==1) adjMat[i,c(i,i+1)] <- c(1,-1)
      else if (i==n) adjMat[i,c(i-1,i)] <- c(-1,1)
      else adjMat[i,c(i-1,i,i+1)] <- c(-1,2,-1)
    }
    return(adjMat)
  }
  # if lags specified, not Griffith style.
  for (i in 1:n) {
    for (l in lags) {
      if (i+l<=n) adjMat[i,(i+l)] <- 1
      if (i-l>0) adjMat[i,(i-l)] <- 1
    }
  }
  return(adjMat)
}

# Function that performs a Moran basis operation given covariates and adjacency matrix
# (see Griffith, 2013; Hughes and Haran, 2013 -- can be accommodated to either setting)
moranEigen_fn <- function(X=NULL, A, k, attractive=T) {
  # PRE: X is n x p covariates (default is NULL -- no orthogonality restriction). 
  # A is the adjacency matrix that the user specifies.
  # POST: Returns a list of k eigenvalues and eigenvectors of the Moran basis matrix.
  if (require(Matrix) && require(RSpectra)) {
    n <- nrow(A)
    if (is.null(X)) X <- rep(1,n) # if NULL, no covariates -- all 1's
    projection <- X %*% solve(t(X) %*% X) %*% t(X) # Projection matrix onto col. space of X
    I <- diag(1,n,n) # identity matrix
    Mx <- t(I-projection) %*% as(A, "dgCMatrix") %*% (I-projection) # Moran expansion on the adj. matrix.
    result <- eigs_sym(k=k, as(Mx, "dgeMatrix")) # Returns a list.
    if (attractive) return(result$vectors[,result$values > 0])
    return(result)
  }
}


#---
# Header for CSV files in dataset
g_header_str <- c("STATION", "STATION_NAME", "ELEVATION", "LATITUDE", "LONGITUDE", "DATE",
  "REPORTTPYE", "HOURLYSKYCONDITIONS", "HOURLYVISIBILITY", "HOURLYPRSENTWEATHERTYPE",
  "HOURLYDRYBULBTEMPF", "HOURLYDRYBULBTEMPC", "HOURLYWETBULBTEMPF", "HOURLYWETBULBTEMPC",
  "HOURLYDewPointTempF", "HOURLYDewPointTempC", "HOURLYRelativeHumidity",
  "HOURLYWindSpeed", "HOURLYWindDirection", "HOURLYWindGustSpeed", "HOURLYStationPressure",
  "HOURLYPressureTendency", "HOURLYPressureChange", "HOURLYSeaLevelPressure",
  "HOURLYPrecip", "HOURLYAltimeterSetting", "DAILYMaximumDryBulbTemp", "DAILYMinimumDryBulbTemp",
  "DAILYAverageDryBulbTemp", "DAILYDeptFromNormalAverageTemp", "DAILYAverageRelativeHumidity",
  "DAILYAverageDewPointTemp", "DAILYAverageWetBulbTemp", "DAILYHeatingDegreeDays",
  "DAILYCoolingDegreeDays", "DAILYSunrise", "DAILYSunset", "DAILYWeather", "DAILYPrecip",
  "DAILYSnowfall", "DAILYSnowDepth", "DAILYAverageStationPressure", "DAILYAverageSeaLevelPressure",
  "DAILYAverageWindSpeed", "DAILYPeakWindSpeed", "PeakWindDirection", "DAILYSustainedWindSpeed",
  "DAILYSustainedWindDirection", "MonthlyMaximumTemp", "MonthlyMinimumTemp",
  "MonthlyMeanTemp", "MonthlyAverageRH", "MonthlyDewpointTemp", "MonthlyWetBulbTemp",
  "MonthlyAvgHeatingDegreeDays", "MonthlyAvgCoolingDegreeDays", "MonthlyStationPressure",
  "MonthlySeaLevelPressure", "MonthlyAverageWindSpeed", "MonthlyTotalSnowfall", "MonthlyDeptFromNormalMaximumTemp",
  "MonthlyDeptFromNormalMinimumTemp", "MonthlyDeptFromNormalAverageTemp", "MonthlyDeptFromNormalPrecip",
  "MonthlyTotalLiquidPrecip", "MonthlyGreatestPrecip", "MonthlyGreatestPrecipDate", "MonthlyGreatestSnowfall",
  "MonthlyGreatestSnowfallDate", "MonthlyGreatestSnowDepth", "MonthlyGreatestSnowDepthDate", "MonthlyDaysWithGT90Temp",
  "MonthlyDaysWithLT32Temp", "MonthlyDaysWithGT32Temp", "MonthlyDaysWithLT0Temp", "MonthlyDaysWithGT001Precip", "MonthlyDaysWithGT010Precip",
  "MonthlyDaysWithGT1Snow", "MonthlyMaxSeaLevelPressureValue", "MonthlyMaxSeaLevelPressureDate", "MonthlyMaxSeaLevelPressureTime",
  "MonthlyMinSeaLevelPressureValue", "MonthlyMinSeaLevelPressureDate", "MonthlyMinSeaLevelPressureTime", "MonthlyTotalHeatingDegreeDays",
  "MonthlyTotalCoolingDegreeDays", "MonthlyDeptFromNormalHeatingDD", "MonthlyDeptFromNormalCoolingDD", "MonthlyTotalSeasonToDateHeatingDD",
  "MonthlyTotalSeasonToDateCoolingDD")

# Reading CSV files with options preset
myRead_csv <- function(file) read.csv(file, header=T, stringsAsFactors=F)

# Reading CSV without header (for merging use)
readNoHeader_csv <- function(file) {
  read.csv(file, header=F, stringsAsFactors=F, skip=1)
}

# Function to merge all CSV files in the path to one data frame
mergeFiles_csv <- function(path) {
  fileNames <- list.files(path, "*.csv", full.names=T)
  allFiles <- lapply(fileNames, readNoHeader_csv)
  result <- do.call(rbind.data.frame, allFiles)
  names(result) <- g_header_str
  return(result)
}

rm(g_header_str, myRead_csv, readNoHeader_csv)
