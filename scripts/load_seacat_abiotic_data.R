 library(tidyverse)
 library(here)
#MEAN SST AND SALINITY AVERAGED ACROSS DEPTH AND MEAN/STAION (MEAN ACROSS HAULS IF  MULTIPLE HAULS)
setwd("~/Documents/GitHub/MAPP/data/ecodaat abiotic data")

#load in all files from this folder that end with csv -- for some reason they  are compiled and have temp but when I  tried  to compile them on my own they didnt work... 

SST_SAL <- list.files(pattern =  "*.csv") %>% 
  purrr::map_df(~read_csv(., col_types = cols("TIME"  = "c"))) %>%
  dplyr::select(1:16) %>% 
  mutate(DEPTH = round(DEPTH)) %>%
  filter(DEPTH < 101) %>% #temp in first 100m, greater  likelihood this impacts ZOOP
  group_by(CRUISE, STATION) %>% #cant use date - there is NAs in many  of the dates. 
  dplyr::summarise(SALINITY = mean(SALINITY), TEMPERATURE = mean(TEMPERATURE)) %>%
  dplyr::rename(#LAT = LATITUDE, LON = LONGITUDE, 
    STATION_NAME = STATION) 



#load in the data that I had to compile from cruises on the server
compiled<-read_csv("compiled from server/compiled_abiotic.csv") %>%
  filter(DEPTH < 101) %>%
  group_by(CRUISE, STATION_NAME) %>%
  dplyr::summarise(SALINITY = mean(SALINITY1), TEMPERATURE = mean(TEMPERATURE1)) 


SST_SAL_complete <- rbind(compiled, SST_SAL)
 
setwd("~/Documents/GitHub/FISH_604")
write_csv(SST_SAL_complete,"data/SST_SAL_complete.csv")
    
    