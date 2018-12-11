# 
# Exploring time series using xts package + visualization
# Data subset: Phoenix airport, AZ, Year 2017.
#

# Required packages
library(tidyverse)
library(GGally)
library(lubridate)
library(xts) # For time series

phoenix_df <- read.csv("Phoenix Airport, AZ.csv", header=T, stringsAsFactors = F)
## date and time are input as strings
parserDate_fn <- function(str) {
  ## PRE: string containing information about date and time (well-formatted)
  ## POST: UTC converted, POSIX vector
  mst <- as.POSIXct(str, tz="MST")
  as_datetime(mst, tz="UTC")
}

phoenix_df <- phoenix_df %>% mutate(DATE = parserDate_fn(DATE))

## Using xts class: subset Year 2017.
# Objective: visualizing diurnal pattern across months.
phWindSpeed_xts <- phoenix_df %>% 
  select(HOURLYWindSpeed) %>% 
  as.xts(order.by = phoenix_df$DATE)

hourlyXtsSummary_fn <- function(xts, subperiod) {
  ## indexes a subperiod (string) of an xts class object
  data <- xts[subperiod]
  temp_df <- data.frame(
    time = index(data),
    windSpeed = as.numeric(data)
    )
  # data.frame of time periods and wind speeds
  summary_df <- temp_df %>% 
    mutate(month = lubridate::month(time),
           hour = lubridate::hour(time)) %>% 
    group_by(month, hour) %>% 
    summarise(avg = mean(windSpeed, na.rm=T),
              std = sd(windSpeed, na.rm=T))
  # Returns data.frame including hourly averages / std's, per each month
  summary_df
}

## Function to visualize patterns for 12 months: use ggplot2
hourlySummaryPlot_fn <- function(df, year, index) {
  # Year for the plot, index is for file name (string)
  plotList <- list()
  for (m in 1:12) {
    m_df <- filter(df, month==m)
    plotList[[m]] <- ggplot(data=m_df, mapping=aes()) + 
      geom_line(mapping=aes(hour, avg), colour="pink") + 
      geom_line(mapping=aes(hour, std), colour="lightblue")
  }
  hourlyPlot_ggm <- ggmatrix(plotList, 2, 6, 
                             yAxisLabels = c("Jan-Jun", "Jul-Dec"),
                             xlab = "Hour", ylab = "Wind Speed",
                             title = paste0("Hourly Average (pink) and Std (blue) across Days in ", year) 
                             ) + theme_classic()
  # Save it as a png file (standard extension for now)
  ggsave(filename=paste0(index, "_wind-hourly-average-std-plot-", year, ".png"), 
         plot=hourlyPlot_ggm, device='png')
}

# Three year, hourly average / std, for each of 12 months
hourlyPh2017 <- hourlyXtsSummary_fn(phWindSpeed_xts, '2017')
hourlyPh2016 <- hourlyXtsSummary_fn(phWindSpeed_xts, '2016')
hourlyPh2015 <- hourlyXtsSummary_fn(phWindSpeed_xts, '2015')

hourlySummaryPlot_fn(hourlyPh2017, 2017, "001.01")
hourlySummaryPlot_fn(hourlyPh2016, 2016, "001.02")
hourlySummaryPlot_fn(hourlyPh2015, 2015, "001.03")
