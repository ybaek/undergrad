# 
# Exploring how to use classes in sp package.
# Data subset: 5 airports in AZ state. Time is of less concern.
#

# Required packages
library(tidyverse)
library(lubridate)
library(GGally)
library(maps)
library(maptools)
library(sp)

# CRS class object for correct projection grid (global)
g_latlongCRS <- CRS("+proj=longlat +datum=WGS84")
# parsing function
parserDate_fn <- function(str) {
  ## PRE: string containing information about date and time (well-formatted)
  ## POST: UTC converted, POSIX vector
  mst <- as.POSIXct(str, tz="MST")
  as_datetime(mst, tz="UTC")
}

#1. Reading spatiotemporal data: 5 airports in AZ state
AZ_df <- read.csv("AZ.csv", header=T, stringsAsFactors=F)
## Data: merged csv file (preprocessed using Bash)
AZ_df <- AZ_df %>% mutate(DATE = parserDate_fn(DATE))
timeGrid_psx <- AZ_df$DATE ## time grid (for spacetime)

#2. map of Arizona counties (space grid)
mapAZ <- maps::map("county", "arizona", plot=F, fill=T)
AZ_IDs <- sapply(strsplit(mapAZ$names, ","), function(x) x[2])
## A character vector of county names
mapAZ_sp <- mapAZ %>% map2SpatialPolygons(IDs = AZ_IDs,
                      proj4string = g_latlongCRS)
## Converting into "SpatialPolygons" class; can produce maps(polygons)

#3. Plotting coordinates of stations on AZ map
## Average wind speed (pretty meaningless; for this purpose)
c_avgWind_df <- AZ_df %>%
  group_by(STATION_NAME) %>% 
  summarise(avg = mean(as.numeric(HOURLYWindSpeed), na.rm=T))

## Prepare a data frame object containing 5 station coordinates.
stationCoords_df <- data.frame(
  station = unique(AZ_df$STATION_NAME),
  long = unique(AZ_df$LONGITUDE),
  lat = unique(AZ_df$LATITUDE),
  avg = c_avgWind_df$avg
)
stationCoords_spdf <- SpatialPointsDataFrame(
  coords = stationCoords_df[2:3],
  data=stationCoords_df[c(1, 4)]
)
## Now converted into "SpatialPointsDataFrame" class
stationCoords_spdf@proj4string <- g_latlongCRS
## @ used for *slots* in S4 objects; proj4string is CRS object

# Overlaying points on a map (base plot)
plot(mapAZ_sp, axes=T, main="Airport locations in AZ")
points(stationCoords_spdf, pch=19)

# Using spplot with scales
spplot(stationCoords_spdf, "avg", 
       scales=list(draw=T)) ## axes to be drawn or no?

# Using ggplot2
ggplot(data=stationCoords_df, 
       mapping=aes(long, lat, colour = avg)) + 
  geom_point() +
  labs(colour="Wind speed (avg)", 
       x="Latitude", y="Longitude",
       title="Airports in AZ") +
  coord_equal()
  