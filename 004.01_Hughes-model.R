#
# Codes implementing spatial mixed model as proposed in
# Hughes and Haran (2013); Hughes (2018)
#

# Preparing data
library(tidyverse)
library(readxl)
library(maps)
library(ngspatial)
source("000_utils.R")
population_tbl <- read_excel("data/PopulationEstimates.xls", sheet=1, skip=2)
unemp_tbl <- read_excel("data/Unemployment.xls", sheet=1, skip=6)
counties_tbl <- read_csv("data/county_adjacency2010.csv", n_max=21721) 
## not reading Puerto Rico/other territories

# Data wrangling to merge tables
# There are in total 3141 county equivalents in U.S. mainland (excluding DC)
DeathMig17_tbl <- population_tbl %>% 
  filter((as.numeric(FIPS) %% 1000)!=0, State!='DC', State!='PR') %>%
  select(State, County = Area_Name, 
         DeathRates = Deaths_2010, NetMig = NET_MIG_2010) %>% 
  mutate(County = paste(County, State, sep=", "))

Unemp17_tbl <- unemp_tbl %>% 
  filter((as.numeric(FIPStxt) %% 1000)!=0, State!='DC', State!='PR') %>%
  select(State, County = Area_name, UnempRates = Unemployed_2017)

uniqueCounties_tbl <- counties_tbl %>% 
  group_by(countyname) %>% 
  summarise(neighbors = paste(neighborname, collapse=";")) %>% 
  ungroup()
  
# Creating an adjacency matrix. Let's focus on a subset for now (Texas: 254 counties).
DeathMigTX_tbl <- DeathMig17_tbl %>% 
  filter(State=="TX")
UnempTX_tbl <- Unemp17_tbl %>% 
  filter(State=="TX")
countiesTX_tbl <- uniqueCounties_tbl %>% 
  filter(grepl('TX', countyname))

c_names <- countiesTX_tbl$countyname
c_n <- nrow(countiesTX_tbl)
adjM <- matrix(0, c_n, c_n)
for (i in 1:c_n) {
  for (j in i:c_n) {
    name <- paste0(";", c_names[j])
    if (grepl(name, countiesTX_tbl$neighbors[i])) adjM[i,j] <- 1
    adjM[j,i] <- adjM[i,j]
  }
  adjM[i,i] <- 0
}

# Response of interest: mortality rate
# Possible covariates: unemployment rate, net migration
dataTX_tbl <- inner_join(DeathMigTX_tbl, UnempTX_tbl, by="County") %>% 
  select(-State.y) %>% 
  rename(State = State.x)

# Since mortality rate is symmetrically distributed, let's fit a gaussian model.
mortalityTX_lm <- lm(DeathRates ~ NetMig + UnempRates, data=dataTX_tbl)
# Using sparse.sglmm
# method="RSR" -> HH2013, "BSF"-> H2018.
mortalityTXHughes_model <- sparse.sglmm(
  DeathRates ~ NetMig + UnempRates, family=gaussian, data=dataTX_tbl,
  A=adjM, method="BSF", tune=list(sigma.s=0.02, sigma.h=0.1), 
  hyper=list(a_h=0.5, b_h=2000)
)




