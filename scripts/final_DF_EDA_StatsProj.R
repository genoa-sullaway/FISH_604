#create full data frame for stats project- then saves data frame, so probbaly dont need to run it again. 

# library(ncdf4)
# library(ncdf4.helpers)
# library(PCICt)
library(here)
library(tidyverse)
library(sf)
library(maps)
library(raster)

##############################################################################################
#load data
full_calanus <- read_csv("data/Calanus_data_14-16.csv")

#get list of unique cruises so you  can bring in temperature data 
#cruise_list <- unique(full_calanus$CRUISE)

##############################################################################################
#load SST and salinity data
#source("scripts/load_seacat_abiotic_data.csv")
SST_SAL_complete<-read_csv("data/SST_SAL_complete.csv")

later<-full_calanus %>%
  filter(YEAR > 2016)

unique(later$CRUISE)

##############################################################################################
#create base map
ak <- sf::st_as_sf(maps::map('world','USA:Alaska',
                             # ylim = c(55,60),
                             #   xlim=c(-140,-130), 
                             plot=FALSE, fill=TRUE))

#create bounds to trim AK map 
bounds_ak <- extent(-176,-150,50,66 )
#bounds_ak <- extent(-164.5,-160, 55,60) 
extent_ak <- st_set_crs(st_as_sf(as(bounds_ak, "SpatialPolygons")), 4326)
ak <- st_intersection(ak, extent_ak) #trim map by intersections 
 plot(ak)
 
#load shape file with ortiz shapes 
ortiz_shp <- st_read("data/ortiz_regions/BSIERP_regions_2012.shp")
ortiz_shp <- subset(ortiz_shp, !DOMAIN == 15)
#Tranform into WGS84 coordinate system
ortiz_shp <- st_transform(ortiz_shp, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

#look at where the domains lie on the map
# ggplot() +
#   geom_sf(data = ortiz_shp, aes(fill = as.character(name)), color = "white")+ #color = "black", fill = "white") +
#   geom_sf(data=ak)

##############################################################################################
#load data and compile into one df
#load climate  data from Litzlow -- Wind etc.  
climate_data <- read_csv("data/climate_data.csv") %>%
  filter(year %in% c(2014, 2015,2016)) %>%
  dplyr::select(year,AO.jfm,SE.wind.Oct.Apr, NW.wind.Oct.Apr, SE.wind.May.Sep, NW.wind.May.Sep,
                summer.cold.pool.extent, ice.area.jfma, north.wind.stress.amj,south.wind.stress.amj) %>%
  dplyr::rename(YEAR = "year")

#combine all covariates
df <- full_calanus %>%
      mutate(DAY = as.numeric(DAY),MONTH = as.numeric(MONTH), YEAR = as.numeric(YEAR) )  %>%
      left_join(SST_SAL_complete) %>% # YAY - THIS WORKS:  c("CRUISE","STATION_NAME")
      left_join(climate_data) 

st_geometry(df) <- NULL
write.csv(df, "data/fulldata_stats_proj.csv") #created this DF on 9/14/2021, saved it so I dont have to keep running this script
##############################################################################################
