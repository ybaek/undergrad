#
# Preprocessing source file
# Creates a sparse adj. matrix between 3141 counties in U.S. mainland.
#

library(tidyverse)
library(Matrix)
library(RSpectra)
# 1. There are in total 3141 counties in the U.S. mainland, excluding DC.
c_counties <- read_csv("data/county_adjacency2010.csv") %>% 
  filter(!str_detect(countyname, "(DC|PR|VI|GU|AS|MP)")) %>% 
  group_by(fipscounty, countyname) %>% 
  summarise(neighbors = paste(fipsneighbor, collapse=";")) %>% 
  ungroup()

## data input errors/up-to-date corrections.
c_counties[87,2] <- "Petersburg Borough, AK"
c_counties[93,2] <- "Kusilvak Census Area, AK"
c_counties[1142,2] <- "LaSalle Parish, LA"
c_counties[1369,2] <- "Otter Tail County, MN"
c_counties[1396,2] <- "Watonwan County, MN"
c_counties[1802,2] <- "Doña Ana County, NM"
c_counties[2417,2] <- "Oglala Lakota County, SD"

c_counties <- c_counties[-2916,] %>% # Bedford City deleted
  arrange(fipscounty, countyname)
g_countyNames <- unique(c_counties$countyname)

# A symmetric binary adjacency matrix
c_n <- nrow(c_counties)
c_codes <- c_counties$fipscounty
A <- matrix(0, c_n, c_n)
for (i in 1:c_n) {
  neighbors <- strsplit(c_counties[i,]$neighbors, ";")[[1]]
  A[i,] <- c_codes %in% neighbors
  A[i,i] <- 0
}
rownames(A) <- unique(c_counties$fipscounty)
colnames(A) <- unique(c_counties$fipscounty)
rm(list=c("c_n", "c_codes", "c_counties"))
