 library(tidyverse)
 library(here)
#MEAN SST AND SALINITY AVERAGED ACROSS DEPTH AND MEAN/STAION (MEAN ACROSS HAULS IF  MULTIPLE HAULS)
#setwd("~/Documents/GitHub/FISH_604")
#load in all files from this folder that end with csv -- for some reason they  are compiled and have temp but when I  tried  to compile them on my own they didnt work... 
temp1<-read_csv("data/ecodaat abiotic data/AE15-01.csv") %>%
  dplyr::rename(#LAT = "GEARINLATITUDE", LON = "GEARINLONGITUDE", 
         STATION_NAME = STATION #, TIME = "GEARIN TIME"
         ) %>%
  dplyr::select(-STATIONID) %>%
  group_by(CRUISE, STATION_NAME) %>%
  dplyr::summarise(SALINITY = mean(SALINITY), TEMPERATURE = mean(TEMPERATURE)) 

temp2<-read_csv("data/ecodaat abiotic data/DY14-08L1_forEcoDAAT.csv") %>%
  dplyr::rename(STATION_NAME = STATION ) %>%
  group_by(CRUISE, STATION_NAME) %>%
  dplyr::summarise(SALINITY = mean(SALINITY), TEMPERATURE = mean(TEMPERATURE)) 


# SST_SAL <- list.files(pattern =  "*.csv") %>% 
#   purrr::map_df(~read_csv(., col_types = cols("TIME"  = "c"))) %>%
#   dplyr::select(1:16) %>% 
#   mutate(DEPTH = round(DEPTH)) %>%
#   filter(DEPTH < 101) %>% #temp in first 100m, greater  likelihood this impacts ZOOP
#   group_by(CRUISE, STATION) %>% #cant use date - there is NAs in many  of the dates. 
#   dplyr::summarise(SALINITY = mean(SALINITY), TEMPERATURE = mean(TEMPERATURE)) %>%
#   dplyr::rename(#LAT = LATITUDE, LON = LONGITUDE, 
#     STATION_NAME = STATION) 

library(here)
#load in the data that I had to compile from cruises on the server
#updated the abiotic data on 11-7-21 to add in 2017-2019 
compiled<-read_csv("data/ecodaat abiotic data/compiled_abiotic.csv") %>%
  filter(DEPTH < 101) %>%
  group_by(CRUISE, STATION_NAME) %>%
  dplyr::summarise(SALINITY = mean(SALINITY1), TEMPERATURE = mean(TEMPERATURE1)) 

# compiled<-read_csv("compiled from server/compiled_abiotic.csv") %>%
#   filter(DEPTH < 101) %>%
#   group_by(CRUISE, STATION_NAME) %>%
#   dplyr::summarise(SALINITY = mean(SALINITY1), TEMPERATURE = mean(TEMPERATURE1)) 

SST_SAL_complete <- rbind(compiled, temp1, temp2) %>%
  dplyr::mutate(CRUISE = case_when(CRUISE == "SQ18-01" ~  "AQ18-01",
         TRUE ~ CRUISE))
 
# temp_list<-data.frame(unique(SST_SAL_complete$CRUISE))

setwd("~/Documents/GitHub/FISH_604")
write_csv(SST_SAL_complete,"data/SST_SAL_complete.csv")
    
    