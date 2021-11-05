#copied from Litzlow github and adapted by GS: https://github.com/mikelitzow/bold-new-pollock
#this script  will provide mean SST and SST Jan-April anomoly for 1x1 (lat and lon degrees) grid. Was going to use it for project but I have station specific data  now. 

library(tidyverse)
library(ncdf4)
library(maps)
library(maptools)
library(mapdata)
library(fields)
library(chron)
library(sf)
library(magrittr)

###################
# add ERSST v5 -- Temperature

download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1950-01-01):1:(2020-2-01)][(0.0):1:(0.0)][(54):1:(65)][(180):1:(205)]", "data/ersst")
#need to make the longitudes in degrees East!!! the negatives mean they are in degrees W. so you add the difference to 180. 

# load and process SST data
nc <- nc_open("data/ersst")

# extract dates

ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))
m <- months(d)
yr <- as.numeric(as.character(years(d)))

# extract study area **** Genoa adapted this from his script!
# 54-60 deg. N, 158-174 deg. E
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1) 

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  #if  I want daily values I tink I change this line

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
# need to drop GOA cells --- dont think I need to do this for the data i selected since i do the spatial  filter later 
# GOA <- c("N54E194", "N54E196", "N54E198", "N54E200", "N54E202", "N56E200", "N56E202")
# SST[,GOA] <- NA

#remove columns that have Nan
#SST <-  SST[ , apply(SST, MARGIN = 2, function(x) sum(is.na(x)) == 0)]

#write_csv(SST, "data/SST_litzlow.csv") #save this df so I dont have to call the URL everytime

####################################### MONTHLY MEAN
year <- c(1950:2019)
month <- c(1:12)
y.m<-expand.grid(year, month) %>%
  arrange(Var1, Var2) %>%
  rename(year = "Var1", month = "Var2")

col = ncol(SST) 
SST.df<-SST %>%
  data.frame() %>%
  tibble::rownames_to_column(var = "date") %>%
  dplyr::select(-date) %>%
  slice(1:(n() - 2)) %>% # remove last two rows
  cbind(y.m) %>%
  gather(1:col, key = "geo_key", value = "temp") %>%
  # gather(1:77, key = "geo_key", value = "temp") %>%
  separate(geo_key, into = c("LAT", "LON"), sep =3) %>%
  separate(LAT, into = c("delete", "LAT"), sep=1) %>%
  separate(LON, into = c("delete", "LON"), sep = 1) %>%
  dplyr::select(year, month, LAT, LON, temp) %>%
  filter(year %in% c(2014,2015,2016)) %>%
  mutate(LON=as.numeric(LON),
         LAT=as.numeric(LAT),
         diff = LON-180,
         LON_new = -1*(180 -diff),
         LON = LON_new) %>% #convert to negative/Eastern degrees
  dplyr::select(-LON_new, -diff)

#plot SST GRID on top of the ortiz shape file 
# sf_SST = st_as_sf(SST.df, coords = c("LON", "LAT"), crs = 4326)
# 
# ggplot() +
#   geom_sf(data = ortiz_shp, aes(fill = as.character(name)), color = "white")+ #color = "black", fill = "white") +
#   geom_sf(data=sf_SST)

#SST.df is now a df with mean SST for that month with lat and long coordinates. but idk how that relates to the ortiz domains. 
#assign SST to spatial regions so it can eventually map onto zoop data
#Load shapefile that Dave sent for Ortiz regions
ortiz_shp <- st_read("data/ortiz_regions/BSIERP_regions_2012.shp")
ortiz_shp <- st_transform(ortiz_shp, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")#Tranform into WGS84 coordinate system
#ortiz_shp <- subset(ortiz_shp, !DOMAIN == 15)
sf_SST_df = st_as_sf(SST.df, coords = c("LON", "LAT"), crs = 4326) # change df to sf
sf_ortiz = st_as_sf(ortiz_shp) #convert your shapefile to sf, too. 

#First do a join that contains the points to filter points that are tootally ouotside the regiono. 
#Then do a NN join to hopefully get a point or each polygon? 
SST_ortiz_df_contains <- st_join(sf_ortiz, sf_SST_df, join = st_contains) 

#need to do some tidying to, rejoin OG lat longs and plot to look at trim
st_geometry(SST_ortiz_df_contains)<-NULL

SST_ortiz_df_contains  %<>%
  dplyr::select(-name, -DOMAIN, -source, -area_km2) %>%
  left_join(SST.df) %>%
  filter(!is.na(LAT))  %>% #one row of NA's
  filter(case_when(LAT == 66  ~ LON > -172, #filtering out these points, but I am not sure why the join didnt work to filter them out on its own, I think it has to do with the weird dateline thing... 
                   LAT == 64 ~ LON > -176,
                   LAT == 62 ~ LON > -178,
                   TRUE ~ LON < -160)) #BASICALLY A FAKE  FILTER TO CLOSE THE CASE WHEN

sf_SST_df_trim = st_as_sf(SST_ortiz_df_contains, coords = c("LON", "LAT"), crs = 4326) # change df to sf

#remove for the plot because of the big lines. 
ortiz_shp_forplot <- subset(ortiz_shp, !DOMAIN == 15)
#plot and check that it trimmed it
ggplot() +
  geom_sf(data = ortiz_shp_forplot, aes(fill = as.character(name)), color = "black", alpha = 0.5) + #color = "black", fill = "white") +
  geom_sf(data=sf_SST_df_trim)

join_with_nn<-st_join(sf_SST_df_trim,sf_ortiz, join = st_nearest_feature) #now that it is trim assign based on nearest feature to make sure each domain gets some data, but the pribolofs still wont get anything.. 

#how did the join perform? 
ggplot() +
  geom_sf(data = ortiz_shp_forplot, aes(fill = as.character(name)), color = "black", alpha = 0.5) + #color = "black", fill = "white") +
  geom_sf(data=join_with_nn, aes(color = name)) 

#summarize to get mean and variance temp per month, year and domain. 
full.MU.SST <- join_with_nn %>%
  dplyr::select(DOMAIN, name, area_km2,year, month,temp) %>%
  group_by(DOMAIN, name, area_km2,year, month) %>%
  dplyr::summarise(mean_temp = mean(temp)) %>%
  rename(YEAR = "year", MONTH = "month")
st_geometry(full.MU.SST) <- NULL
######################################################################################
#calculate anomolies

#this is used to calculate monthly means across the time series so that it can be compared to each year/month and calculate
#the anomoly for that month compared to the average on that whole time series! 
f <- function(x) tapply(x, m, mean)
#I am not sure why this part of the code says 1981 when it starts at 1950?? check this!
# mu.SST <- t(as.matrix(colMeans(SST[yr %in% 1981:2010,])))
# rep(c(1:12),times=length(d))

mu.SST <- apply(SST, 2, f)	# Compute monthly means for each time series (location)
#mu.SST <- apply(SST[yr %in% 1981:2010,], 2, f)	# Compute monthly means for each time series (location)
mu.SST <- mu.SST[rep(1:12, floor(length(d)/12)),] 
 
SST<- SST[-(841:842),]

anom.SST <- SST - mu.SST #get monthly anomolies

col = ncol(anom.SST) 
anom.SST.df<-anom.SST %>%
  data.frame() %>%
  cbind(y.m) %>%
  gather(1:col, key = "geo_key", value = "anom_temp") %>%
  filter(month <5) %>% # filter to just have Jan - April anomolies
  separate(geo_key, into = c("LAT", "LON"), sep =3) %>%
  separate(LAT, into = c("delete", "LAT"), sep=1) %>%
  separate(LON, into = c("delete", "LON"), sep = 1) %>%
  dplyr::select(year, month, LAT, LON, anom_temp) %>%
  filter(year %in% c(2014,2015,2016)) %>%
  mutate(LON=as.numeric(LON),
         LAT=as.numeric(LAT),
         diff = LON-180,
         LON_new = -1*(180 -diff),
         LON = LON_new) %>% #convert to negative Eastern degrees
  dplyr::select(-LON_new, -diff) 

#filter with shape files -- [this filters out a lot of the NAs] and then get mean Jan - April anomolies.  
#change to sf
sf_anom_SST_df = st_as_sf(anom.SST.df, coords = c("LON", "LAT"), crs = 4326) # change df to sf
#do the first join
anom_SST_ortiz_df_contains <- st_join(sf_ortiz,sf_anom_SST_df, join = st_contains)
st_geometry(anom_SST_ortiz_df_contains)<-NULL

#first trim is done - but still have to take out a few points, add lat long back in for the NN trim. 
#test <- anom_SST_ortiz_df_contains %>%
anom_SST_ortiz_df_contains  %<>%
  dplyr::select(-name, -DOMAIN, -source, -area_km2) %>%
  left_join(anom.SST.df) %>% #add lat longs back in
  filter(!is.na(LAT))  %>% #one row of NA's
  filter(case_when(LAT == 66  ~ LON > -172, #filtering out these points, but I am not sure why the join didnt work to filter them out on its own, I think it has to do with the weird dateline thing... 
                   LAT == 64 ~ LON > -176,
                   LAT == 62 ~ LON > -178,
                   TRUE ~ LON < -160)) #BASICALLY A FAKE  FILTER TO CLOSE THE CASE WHEN
#make DF into a SF
sf_anom_SST_df = st_as_sf(anom_SST_ortiz_df_contains, coords = c("LON", "LAT"), crs = 4326) # change df to sf

#nearest neighbor join
join_anom_nn<-st_join(sf_anom_SST_df,sf_ortiz, join = st_nearest_feature) #now that it is trim assign based on nearest feature to make sure each domain gets some data, but the pribolofs still wont get anything.. 

#plot and check that it trimmed it
ggplot() +
  geom_sf(data = ortiz_shp_forplot, aes(fill = as.character(name)), color = "black", alpha = 0.5) + #color = "black", fill = "white") +
  geom_sf(data=join_anom_nn, aes(color = name)) 

#sf_SST_df_trim = st_as_sf(sf_anom_SST_df, coords = c("LON", "LAT"), crs = 4326) # change df to sf

full.ANOM.SST <- join_anom_nn  %>%
    group_by(DOMAIN, name, area_km2,year) %>%
    summarise(mean_temp_anom_JanApril = mean(anom_temp)) %>% #mean jan-april anomoly per domain!
    rename(YEAR = "year")
st_geometry(full.ANOM.SST) <- NULL

#I dont think I need the rest of this because it is already in climate data
#############################################################################
#plot data spatially  - some temp variability early and late in the year but in summer it looks hot every where
ggplot(ak) +
  geom_sf() +
  # coord_sf(#crs=4326, #transforms it 
  #          #xlim = mapRange[c(1:2)], 
  #          #ylim = mapRange[c(3:4)], # notice what happens when y
  #          # ylim = c(40,45),
  #          # xlim=c(100,110),
  #          expand = TRUE,
  #          clip = "off")+
  geom_point(data = SST.df, aes(x=LON, y=LAT, color =temp))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(year~month)+
  scale_colour_gradient(low = "blue", high = "red")
#############################################################################
#look at SST.df and make sure you understand all the temporal coverage. There are some blanks in the plotted data that I dont think should be there...
#YUP THIS DATA PROVIDE FULL TEMPORAL COVERAGE FOR 12 MONTHS AND 3 YEARS 
#how many months have temp data within each year and station
month_temporal_check <- SST.df %>%
  group_by(year, LAT, LON) %>%
  count(month) %>%
  group_by(year, LAT, LON) %>%
  summarise(monthsum = sum(n)) #there are 12 months within each year and lat lon... 

year_temporal_check <- SST.df %>%
  group_by(LAT, LON) %>%
  count(year) %>%
  group_by(LAT, LON) %>%
  summarise(yearsum = sum(n)) #should be 12x3 for all n's
  



########################################################################
# # add AO
# 
# ao <- read.csv("data/ao.csv")
# head(ao)
# 
# # restrict to JFM and get annual means
# ao <- ao %>%
#   filter(month <= 3)
# 
# ao <- tapply(ao$value, ao$year, mean)
# 
# # limit to 1951-2019
# ao <- ao[names(ao) %in% 1951:2019]
# 
# clim.dat$AO.jfm <- ao
# 
# ########
# # now NCEP/NCAR winds
# # first U-wind (zonal / east-west winds)
# URL <- 
#   "http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_e77d_1b03_9908.nc?uwnd[(1948-01-01):1:(2020-01-01T00:00:00Z)][(45):1:(75)][(150):1:(225)]"
# 
# download.file(URL, "data/NCEP.NCAR.u-wind.nc")
# 
# # and test
# test <- nc_open("data/NCEP.NCAR.u-wind.nc")
# test
# 
# x <- ncvar_get(test, "longitude")
# y <- ncvar_get(test, "latitude")
# uwnd <- ncvar_get(test, "uwnd", verbose = F)
# 
# # Change data into a matrix with months / cells for rows / columns
# uwnd <- aperm(uwnd, 3:1)  
# uwnd <- matrix(uwnd, nrow=dim(uwnd)[1], ncol=prod(dim(uwnd)[2:3]))  
# 
# z <- colMeans(uwnd, na.rm=T) # mean value for each cell
# z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# image(x,y,z, col=tim.colors(64), xlab = "", ylab = "")
# 
# contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
# map('world2Hires', add=T, lwd=1)
# # looks good!
# 
# # first V-wind (meridional / north-south winds)
# URL <- 
#   "http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_aa14_316a_e154.nc?vwnd[(1948-01-01):1:(2020-01-01T00:00:00Z)][(45):1:(75)][(150):1:(225)]"
# 
# download.file(URL, "data/NCEP.NCAR.v-wind.nc")
# 
# # and test
# test <- nc_open("data/NCEP.NCAR.v-wind.nc")
# test
# 
# x <- ncvar_get(test, "longitude")
# y <- ncvar_get(test, "latitude")
# vwnd <- ncvar_get(test, "vwnd", verbose = F)
# 
# # Change data into a matrix with months / cells for rows / columns
# vwnd <- aperm(vwnd, 3:1)  
# vwnd <- matrix(vwnd, nrow=dim(vwnd)[1], ncol=prod(dim(vwnd)[2:3]))  
# 
# z <- colMeans(vwnd, na.rm=T) # mean value for each cell
# z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
# 
# contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
# map('world2Hires', add=T, lwd=1)
# # looks good!
# 
# #####################
# # now we need daily winds for 60ºN 170ºW to calculate proportion of NW/SE winds
# 
# URL <- "http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_a66c_3524_57cf.nc?uwnd[(1948-01-01):1:(2019-12-31T00:00:00Z)][(60):1:(60)][(190):1:(190)]"
# 
# download.file(URL, "data/NCEP.NCAR.daily.u-wind.nc")
# 
# # and load
# dat <- nc_open("data/NCEP.NCAR.daily.u-wind.nc")
# 
# uwnd <- ncvar_get(dat, "uwnd", verbose = F)
# 
# # extract dates
# raw <- ncvar_get(dat, "time") # seconds since 1-1-1970
# h <- raw/(24*60*60)
# d <- dates(h, origin = c(1,1,1970))
# m <- months(d)
# yr <- as.numeric(as.character(years(d)))
# 
# # add v-wind
# URL <- "http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_3fcd_f037_d8fc.nc?vwnd[(1948-01-01):1:(2019-12-31T00:00:00Z)][(60):1:(60)][(190):1:(190)]"
# 
# download.file(URL, "data/NCEP.NCAR.daily.v-wind.nc")
# dat <- nc_open("data/NCEP.NCAR.daily.v-wind.nc")
# vwnd <- ncvar_get(dat, "vwnd", verbose = F)
# 
# # check that the dates are identical between the two data sets
# raw <- ncvar_get(dat, "time") # seconds since 1-1-1970
# h <- raw/(24*60*60)
# d.v <- dates(h, origin = c(1,1,1970))
# 
# identical(d, d.v) # True!
# 
# # now make a dataframe with date info and uwnd
# 
# daily.wind <- data.frame(date=as.character(d),
#                          month=m,
#                          year=yr,
#                          uwnd=uwnd,
#                          vwnd=vwnd)
# 
# # calculate daily direction from u- and v- vectors
# daily.wind$direction <-(180/pi * atan2(-daily.wind$uwnd, -daily.wind$vwnd))+180
# range(daily.wind$direction) # perfecto
# 
# # now make columns of 1s and 0s to indicate if the wind is blowing NW or SE
# 
# daily.wind$NW <- daily.wind$SE <- 0
# 
# for(i in 1:nrow(daily.wind)){
#   # i <- 2
#   if(daily.wind$dir[i] >=105 & daily.wind$dir[i] <=165) daily.wind$SE[i] <- 1
#   if(daily.wind$dir[i] >=285 & daily.wind$dir[i] <=345) daily.wind$NW[i] <- 1
# }
# 
# # function to calculate the proportion of days with a particular wind direction
# f <- function(x) sum(x)/sum(!is.na(x)) 
# 
# # get monthly sums of proportion of days with wind from each direction
# prop.SE <- tapply(daily.wind$SE, list(yr,m), f)
# prop.NW <- tapply(daily.wind$NW, list(yr,m), f)
# 
# # Danielson et al. 2012 GRL recommend seasons of Oct-Apr and May-Sept
# 
# MaySepNW <- prop.NW[,5:9] 
# MaySepSE <- prop.SE[,5:9]
# 
# # now get summer means
# sumNW <- rowMeans(MaySepNW)
# plot(names(sumNW), sumNW, type="o")
# sumSE <- rowMeans(MaySepSE)
# plot(names(sumSE), sumSE, type="o")
# 
# # and winter...
# OctAprNW <- prop.NW[2:nrow(prop.NW), 1:4]
# #add in Oct Nov Dec from previous year
# OctAprNW <- cbind(OctAprNW, prop.NW[1:(nrow(prop.NW)-1),10:12]) 
# colMeans(OctAprNW, na.rm=T) # check how the months compare - pretty similar
# 
# # get means
# winNW <- rowMeans(OctAprNW)
# plot(names(winNW), winNW, type="o") # more of a coherent trend than the summer TS!
# 
# OctAprSE <- prop.SE[2:nrow(prop.SE), 1:4]
# #add in Oct Nov Dec from previous year
# OctAprSE <- cbind(OctAprSE, prop.SE[1:(nrow(prop.SE)-1),10:12]) 
# colMeans(OctAprSE, na.rm=T) # check how the months compare - pretty similar
# 
# # get means
# winSE <- rowMeans(OctAprSE)
# plot(names(winSE), winSE, type="o") 
# 
# # very unusual proportions of winter proportions in recent years! should double-check these!
# 
# # add to climate data
# 
# clim.dat$SE.wind.Oct.Apr <- winSE[names(winSE) %in% clim.dat$year]
# clim.dat$NW.wind.Oct.Apr <- winNW[names(winNW) %in% clim.dat$year]
# 
# # and summer
# 
# clim.dat$SE.wind.May.Sep <- sumSE[names(sumSE) %in% clim.dat$year]
# clim.dat$NW.wind.May.Sep <- sumNW[names(sumNW) %in% clim.dat$year]
# 
# #############
# # now add summer bottom trawl temperatures and cold pool extent
# 
# dat <- read.csv("data/annual.environmental.data.csv")
# 
# head(dat) # note that these are only from 1988, when more stations in the northern EBS were added - could potentially add earlier years...
# # also note that we are not accounting for diferences in the seasonal timing of sampling in different years!!
# 
# clim.dat$summer.bottom.temp <- clim.dat$summer.cold.pool.extent <- NA
# 
# clim.dat$summer.bottom.temp <- dat$AVG_BT[match(clim.dat$year, dat$Year)]
# clim.dat$summer.cold.pool.extent <- dat$CP_EXTENT[match(clim.dat$year, dat$Year)]
# 
# ###########################
# # and sea ice
# # will consider two different kinds of data here - 
# # March ice concentration at PMEL moorings data
# # and also...monthly Bering Sea ice-covered area anomalies from NSIDC: https://nsidc.org/data/NSIDC-0192/versions/3
# 
# # open the area data (non-anomalies) and examine
# dat <- read.csv("data/gsfc.bootstrap.month.area.1978-2018.n.csv")
# 
# head(dat)
# dat$dec.yr <- dat$Year+(dat$Mon-0.5)/12
# 
# ggplot(dat, aes(dec.yr, Bering)) +
#   geom_line()
# 
# # pretty wild!
# 
# # now look at climatology
# mean.dat <- dat %>%
#   group_by(Mon) %>%
#   summarise(mean=mean(Bering))
# 
# ggplot(mean.dat, aes(Mon, mean)) +
#   geom_line() +
#   geom_point()
# 
# # so I think we'll average anomalies over JFMA - that's the highest extent period
# 
# # load area data 
# dat <- read.csv("data/gsfc.bootstrap.month.anomaly.area.1978-2018.n.csv")
# 
# head(dat)
# 
# 
# # restrict to JFMA
# dat <- dat %>%
#   filter(Mon %in% 1:4)
# 
# # make sure there are no 0 ice months as that would read as average conditions!
# sum(dat$Bering==0) #a-ok!
# 
# total.ice <- tapply(dat$Bering, dat$Year, mean)
# 
# # plot to check
# plot(names(total.ice), total.ice, type="l")
# 
# # add to clim.dat
# clim.dat$ice.area.jfma <- NA
# clim.dat$ice.area.jfma <- total.ice[match(clim.dat$year, names(total.ice))]
# 
# # now add March % ice concentration at moorings sites
# dat <- read.csv("data/March mooring site ice concentration.csv")
# 
# head(dat)
# names(dat)[2:4] <- c("M4", "M5", "M8")
# 
# clim.dat$m4.march.ice <- clim.dat$m5.march.ice <- clim.dat$m8.march.ice <- NA
# 
# clim.dat$m4.march.ice <- dat$M4[match(clim.dat$year, dat$Date)]
# clim.dat$m5.march.ice <- dat$M5[match(clim.dat$year, dat$Date)]
# clim.dat$m8.march.ice <- dat$M8[match(clim.dat$year, dat$Date)]
# 
# # and maisie Bering Sea ice cover - it's a repeat of the ice area time series,
# # and only begins in 2006, so not as good, but is also updated in near-real time, 
# # so gives us more recent information...
# 
# # downloaded from https://nsidc.org/data/masie/
# 
# dat <- read.csv("data/maisie.csv")
# 
# head(dat)
# dat$year <- floor(dat$yyyyddd/1000)
# dat$day <- 1000*(dat$yyyyddd/1000 - dat$year)
# 
# # limit to day 1-120
# # and drop 2020, which isn't yet complete
# 
# dat <- dat %>% 
#   filter(day <= 120, year <= 2019)
# 
# names(dat)
# 
# dat <- dat[,c(14,19)]
# 
# maisie <- tapply(dat[,1], dat$year, mean)
# 
# clim.dat$maisie.ice.extent.jfma <- NA
# 
# clim.dat$maisie.ice.extent.jfma <- maisie[match(clim.dat$year, names(maisie))]
# 
# ########################
# # now! NCEP/NCAR wind stress
# 
# # identify latest year and month needed
# year <- 2019
# month <- "12"
# query <- c("e025_be03_a4bb.nc?vflx",
#            "803d_41d9_b553.nc?uflx")
# 
# variable <- c("vflx", "uflx")
# 
# for(i in 1:length(query)){
#   URL <- paste("http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_", query[i], "[(1948-01-01):1:(", year, "-",
#                month, "-01T00:00:00Z)][(19.99970054626):1:(69.52169799805)][(120):1:(249.375)]", sep="")
#   
#   download.file(URL, paste("data/North.Pacific.NCEP.NCAR.", variable[i], sep=""))
# }
# 
# # and load/process
# dat <- nc_open("data/North.Pacific.NCEP.NCAR.uflx")
# dat
# 
# x <- ncvar_get(dat, "longitude")
# y <- ncvar_get(dat, "latitude")
# uflx <- ncvar_get(dat, "uflx", verbose = F)
# 
# # extract dates
# raw <- ncvar_get(dat, "time") # seconds since 1-1-1970
# h <- raw/(24*60*60)
# d <- dates(h, origin = c(1,1,1970))
# m <- months(d)
# yr <- as.numeric(as.character(years(d)))
# 
# # change to matrix
# uflx <- aperm(uflx, 3:1)  
# uflx <- matrix(uflx, nrow=dim(uflx)[1], ncol=prod(dim(uflx)[2:3]))  
# 
# # make vectors of lat/long and add (with date) as dimnames
# lat <- rep(y, length(x))   
# lon <- rep(x, each = length(y)) 
# 
# dimnames(uflx) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
# 
# # now limit to area of interest!!
# keep <- y > 54 & y < 67
# y <- y[keep]
# 
# keep <- x > 187 & x < 202
# x <- x[keep]
# 
# keep.lon <- lon %in% x
# 
# keep.lat <- lat %in% y
# uflx <- uflx[,keep.lat & keep.lon]
# 
# 
# # now remove land and other cells we don't want
# uflx[,c(1,2,3,7,8,28,29,33,35,36,39,40,42,43,46:51,53:56)] <- NA 
# 
# 
# z <- colMeans(uflx, na.rm=T)
# z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# image(x,y,z, col=tim.colors(64), xlab = "", ylab = "")
# 
# contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
# map('world2Hires', add=T, lwd=1)
# 
# # looks good - same routine for v-stress!
# dat <- nc_open("data/North.Pacific.NCEP.NCAR.vflx")
# dat
# 
# x <- ncvar_get(dat, "longitude")
# y <- ncvar_get(dat, "latitude")
# vflx <- ncvar_get(dat, "vflx", verbose = F)
# 
# # change to matrix
# vflx <- aperm(vflx, 3:1)  
# vflx <- matrix(vflx, nrow=dim(vflx)[1], ncol=prod(dim(vflx)[2:3]))  
# 
# # make vectors of lat/long and add (with date) as dimnames
# lat <- rep(y, length(x))   
# lon <- rep(x, each = length(y)) 
# 
# dimnames(vflx) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
# 
# # now limit to area of interest!!
# keep <- y > 54 & y < 67
# y <- y[keep]
# 
# keep <- x > 187 & x < 202
# x <- x[keep]
# 
# keep.lon <- lon %in% x
# 
# 
# keep.lat <- lat %in% y
# lat <- lat[keep.lat & keep.lon]
# vflx <- vflx[,keep.lat & keep.lon]
# 
# # now remove land and other cells we don't want
# vflx[,c(1,2,3,7,8,28,29,33,35,36,39,40,42,43,46:51,53:56)] <- NA 
# 
# z <- colMeans(vflx, na.rm=T)
# z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# image(x,y,z, col=tim.colors(64), xlab = "", ylab = "")
# 
# contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
# map('world2Hires', add=T, lwd=1)
# 
# # get total stress
# total.stress <- sqrt(uflx^2 + vflx^2)
# 
# # limit to spring/summer - AMJ
# keep <- m %in% c("Apr", "May", "Jun")
# yr <- yr[keep]
# 
# total.stress <- total.stress[keep,]
# 
# # limit to north/south
# north.stress <- rowMeans(total.stress[,lat>60], na.rm=T)
# south.stress <- rowMeans(total.stress[,lat<60], na.rm=T)
# 
# north.stress <- tapply(north.stress, yr, mean)
# plot(names(north.stress), north.stress, type="l")
# 
# south.stress <- tapply(south.stress, yr, mean)
# plot(names(south.stress), south.stress, type="l")
# 
# clim.dat$north.wind.stress.amj <- north.stress[match(clim.dat$year, names(north.stress))]
# clim.dat$south.wind.stress.amj <- south.stress[match(clim.dat$year, names(south.stress))]
# 
# ########################
# # plot clim.dat to check
# 
# plot.dat <- clim.dat %>%
#   pivot_longer(-year, names_to = "key", values_to = "value")
# 
# ggplot(plot.dat, aes(year, value)) +
#   geom_line() +
#   facet_wrap(~key, scales="free_y")
# 
# ggsave("figs/climate time series collected to date.png", width=10, height=8, units='in')
# 
# # save!
# 
# write.csv(clim.dat, "data/climate data.csv", row.names = F)