library(here)
library(tidyverse)
library(sf)
library(maps)
library(raster)

#copied form MAPP foder but will use it for the stats project

# Goal is to make a DF with all species of interest that has been correctly filtered for a time series
# I used "figure_out_zoop_filter.R" to understand what the unique combinations were for each species/through time etc. 
# Used Calanus_Zoop.RmD from DK to get an initial understanding of what needs to be tidyed in the data

####  NOTES:
# Currently the multinet and tucker trawl are filtered out because these survey at specific depths and I want depth integrated. 
# Currently keep all life stages in the data frame since ROMS has aggregated across life stages that include 1 and 2, but review this with DK because methods for early lifestages may be squirelly. 
# No c1's in the early years for Neocalanus and Calanus?? -- is this a methods thing? 
################################################################################################################
#Connect to oracle database and download the latest file, save to computer. you may not need to do this every time
################################################################################################################

# #Create connect to the AFSC database. This opens the channel to the Oracle database and provide pasword userID and password. This needs to be set up in advance to work.
# 
# AFSC_Connect <- odbcConnect("AFSC", uid="KIMMELD", pwd="4Scottolana3!"
#                             , believeNRows=FALSE)
# 
# #Delete and refresh table to draw from, in this case it is SPECIMEN_MAIN_GEOM. The table is deleted and rebuilt using this code, that way any updates pushed to the database are captured here by rebuilding the table each time
# 
# sqlQuery(AFSC_Connect,"DROP TABLE SPECIMEN_MAIN_GEOM;")
# 
# sqlQuery(AFSC_Connect,"CREATE TABLE SPECIMEN_MAIN_GEOM AS SELECT * FROM ECODAAT.SPECIMEN_MAIN_GEOM;")
# 
# #Query the database. This section runs the SQL query to build the dataset and add it to the Global Environment. This is a SQL script. One could alter this script to pull data in different ways from the database. Here, I'm just pulling all of the data for manipulation in R
# 
# zoopdata <- sqlQuery(AFSC_Connect, "SELECT BOTTOM_DEPTH,
# CRUISE, DAY, DIS_PERVOLM2, DIS_PERVOLM3, EST_NUM_PERM2, EST_NUM_PERM3, FOCI_ID, FOCI_SAMPLE_ID, GEAR_NAME, 
# GEOGRAPHIC_AREA, GMT_DATE_TIME_TXT, HAUL_ID, HAUL_NAME, HAUL_PERFORMANCE, LAT, LON, MAX_GEAR_DEPTH, MESH,
# MIN_GEAR_DEPTH, MONTH, NET, SAMPLE_DEPTH, SEX, SEX_NAME, SIZE_NAME, SPECIMEN_FORM, STAGE, STAGE_NAME, STATION_NAME,
# TAXON_NAME, TAXON_SIZE, VOLUME_FILTERED, YEAR, ZOOP_COPEPOD_NAUPLII, ZOOP_EUPHAUSIID_EGG
# FROM SPECIMEN_MAIN_GEOM WHERE ORIG_DB LIKE 'BOB';", stringsAsFactors=FALSE)
# 
# #Close database connection so you don't remain connected to the server
# 
# odbcClose(AFSC_Connect)
# 
# 
################################################################################################################
################################################################################################################
# 
##Recode mesh sizes for 150, 333, 500. There are some errors in mesh size in the database, so I just correct for some of these mistakes and to make the mesh sizes similar (sometimes 150 or 500 is entered, when they should be standardized as 153 adn 505)
# 
# zoopdata$MESH[zoopdata$MESH==150] <- 153
# zoopdata$MESH[zoopdata$MESH==154] <- 153
# zoopdata$MESH[zoopdata$MESH==1153] <- 153
# zoopdata$MESH[zoopdata$MESH==335] <- 333
# zoopdata$MESH[zoopdata$MESH==500] <- 505
# 
# write.csv(zoopdata, "AllZoopRaw.csv")

################################################################################################################
#start here if you are just using data from computer
################################################################################################################
sp_list <- c("Calanus marshallae", "Calanus pacificus", "Calanus glacialis",
"Neocalanus spp.",
"Pseudocalanus spp.",
"Thysanoessa raschii",
"Thysanoessa inermis")

zoopdata <- read.csv("data/AllZoopRaw.csv")  %>%
  #filter(!GEAR_NAME %in% c("TUCK1", "MULTINET")) %>%
  mutate(TAXON_NAME = case_when(TAXON_NAME %in% c("Neocalanus cristatus", "Neocalanus flemingeri", "Neocalanus flemingeri / plumchrus", 
                                                  "Neocalanus plumchrus") ~ "Neocalanus spp.", # N flemingeri and N plumchrus were not speciated in old methods
                                TAXON_NAME %in% c("Pseudocalanus spp. CI-CIII") ~ "Pseudocalanus spp.",
                                TRUE ~ TAXON_NAME),
         GEAR_NAME =case_when(GEAR_NAME == "V60BON" ~ "60BON", #not a  typo, this is a gear type. estimates of abundance should be the same. 
                              
                              TRUE~GEAR_NAME)) %>%
  filter(!HAUL_PERFORMANCE == "FAIL",
         GEAR_NAME %in% c("20BON", "60BON")) %>% #pull out just the bongo net gear type. 
filter(!GEOGRAPHIC_AREA == "GOA")
#unique(zoopdata$GEOGRAPHIC_AREA)

#gear notes
#LGCB -  small net inside Tucker, == 20Bon, 
#Sled - epibenthic sampling, cant integrate over depth. 
#calvet - 53 micron, vertical tow, microzoop protocol, good for naupli
#quad net - 20bon, 4 nets. could be useful.  
#ctdb - CTD bottle, microzoop
#methot - 3x3 square, for euphasiids. probably not many, keep for Euphasiids. 
#multinet - mocness, 1 tunnel that opens and closes nets. depth specific
# ikmt - fish trawl, midwater. fine mesh. 
unique(zoopdata$GEAR_NAME)

################################################################################################################
#filter data based on #figure_out_zoop_filter.R and zoop_methods_table.xlsx to create a timeseries for species of interest
################################################################################################################
#unique(zoop_species_interest$STAGE_NAME) #easy to copy and paste into the code below 

neocalanus_spp <- c("Neocalanus spp.", "Neocalanus flemingeri", "Neocalanus flemingeri / plumchrus","Neocalanus plumchrus")
pseudocalanus <- c("Pseudocalanus spp.","Pseudocalanus spp. CI-CIII")
euphasiids<- c("Thysanoessa raschii", "Thysanoessa inermis")

zoop_filtered<- zoopdata %>% 
  filter(case_when(
## CALANUS  MARSHALLAE
    #Old 80's - 2010
    TAXON_NAME=="Calanus marshallae" & YEAR < 2011 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)") ~ SPECIMEN_FORM == "C",
    TAXON_NAME=="Calanus marshallae" & YEAR < 2011 & STAGE_NAME %in% c("C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                        "A + J (ADULT/JUVENILE)","JUVENILE", "ADULT" ) ~ SPECIMEN_FORM == "B",
    #Middle 2011-2018
    TAXON_NAME=="Calanus marshallae" & YEAR > 2010 & YEAR < 2019 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)") ~ SPECIMEN_FORM %in% c("C", "H", "G"), 
    TAXON_NAME=="Calanus marshallae" & YEAR > 2010 & YEAR < 2019 & STAGE_NAME %in% c("C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                    "A + J (ADULT/JUVENILE)","JUVENILE", "ADULT") ~ SPECIMEN_FORM %in% c("B","G"), 

    #Recent 2019 - present
    TAXON_NAME=="Calanus marshallae" & YEAR > 2018 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)") ~ SPECIMEN_FORM %in% c("L", "H"),
    TAXON_NAME=="Calanus marshallae" & YEAR > 2018 & STAGE_NAME %in% c("C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                                  "A + J (ADULT/JUVENILE)","JUVENILE","ADULT") ~ SPECIMEN_FORM %in% c("K","L"),
## CALANUS  PACIFICUS - relatively rare in bering. question if they become more prevalent during warming.
    #Old 80's - 2010  
        #no stage 1-3 in the old dataset
    TAXON_NAME=="Calanus pacificus" & YEAR < 2011 & STAGE_NAME %in% c("C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                                                       "A + J (ADULT/JUVENILE)","JUVENILE", "ADULT" ) ~ SPECIMEN_FORM == "B",
    #Middle 2011-2018
    TAXON_NAME=="Calanus pacificus" & YEAR > 2010 & YEAR < 2019 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)") ~ SPECIMEN_FORM %in% c("H"), 
    TAXON_NAME=="Calanus pacificus" & YEAR > 2010 & YEAR < 2019 & STAGE_NAME %in% c("C - 3 (COPEPODITE III)")  ~  SPECIMEN_FORM %in% c("G"), 
    TAXON_NAME=="Calanus pacificus" & YEAR > 2010 & YEAR < 2019 & STAGE_NAME %in% c("C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                                                                     "A + J (ADULT/JUVENILE)","JUVENILE", "ADULT") ~ SPECIMEN_FORM %in% c("B" ), 
    #Recent 2019 - present
    TAXON_NAME=="Calanus pacificus" & YEAR > 2018 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)") ~ SPECIMEN_FORM %in% c("H"),
    TAXON_NAME=="Calanus pacificus" & YEAR > 2018 & STAGE_NAME %in% c("C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                                                       "A + J (ADULT/JUVENILE)","JUVENILE","ADULT") ~ SPECIMEN_FORM %in% c("K"),
## CALANUS  GLACIALIS # added for arctic work. often combined to calanus spp. 
    #Old 80's - 2010  ---- NO GLACIALIS ENUMERATED BEFORE 2011.  
    #Middle 2011-2018
    TAXON_NAME=="Calanus glacialis" & YEAR > 2010 & YEAR < 2019 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)") ~ SPECIMEN_FORM %in% c("H"), 
    TAXON_NAME=="Calanus glacialis" & YEAR > 2010 & YEAR < 2019 & STAGE_NAME %in% c("C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                                                                    "A + J (ADULT/JUVENILE)","JUVENILE", "ADULT")  ~  SPECIMEN_FORM %in% c("G"), 
    #Recent 2019 - present
    TAXON_NAME=="Calanus glacialis" & YEAR > 2018 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)") ~ SPECIMEN_FORM %in% c("L", "H"),
    TAXON_NAME=="Calanus glacialis" & YEAR > 2018 & STAGE_NAME %in% c("C - 3 (COPEPODITE III)") ~ SPECIMEN_FORM %in% c("G", "L"),  
    TAXON_NAME=="Calanus glacialis" & YEAR > 2018 & STAGE_NAME %in% c("C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)",
                                                                      "A + J (ADULT/JUVENILE)","JUVENILE","ADULT") ~ SPECIMEN_FORM %in% c("K"),

## NEOCALANUS #cristatus should be stand alone filter, cristatus is very large. keep it separated and see if it fits into ROMS. 
    #Old 80's - 2010
    TAXON_NAME %in% neocalanus_spp & YEAR < 2011 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)") ~ SPECIMEN_FORM == "C",
    TAXON_NAME %in% neocalanus_spp & YEAR < 2011 & STAGE_NAME %in% c("C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                                                       "A + J (ADULT/JUVENILE)","JUVENILE", "ADULT" ) ~ SPECIMEN_FORM == "B",
    #Middle 2011-2018
    TAXON_NAME %in% neocalanus_spp & YEAR > 2010 & YEAR < 2019 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)") ~ SPECIMEN_FORM %in% c("H", "G"), 
    TAXON_NAME %in% neocalanus_spp & YEAR > 2010 & YEAR < 2019 & STAGE_NAME %in% c("C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                                                                     "A + J (ADULT/JUVENILE)","JUVENILE", "ADULT") ~ SPECIMEN_FORM %in% c("B"), 
    #Recent 2019 - present
    TAXON_NAME %in% neocalanus_spp & YEAR > 2018 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)") ~ SPECIMEN_FORM %in% c("L", "H"),
    TAXON_NAME %in% neocalanus_spp & YEAR > 2018 & STAGE_NAME %in% c("C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                                                       "A + J (ADULT/JUVENILE)","JUVENILE","ADULT") ~ SPECIMEN_FORM %in% c("K"),
    
## PSEUDOCALANUS # check on the 60BON for early year-- see what was happening here. 
    #Old 80's - 2010
    TAXON_NAME %in% pseudocalanus & YEAR < 2011 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)","C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                                                    "A + J (ADULT/JUVENILE)","JUVENILE", "ADULT" ) ~ SPECIMEN_FORM == "C",
    #Middle 2011-2018
    TAXON_NAME %in% pseudocalanus & YEAR > 2010 & YEAR < 2019 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)","C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                                                                  "A + J (ADULT/JUVENILE)","JUVENILE", "ADULT") ~ SPECIMEN_FORM %in% c("H", "C"), 
    #Recent 2019 - present
    TAXON_NAME %in% pseudocalanus & YEAR > 2018 & STAGE_NAME %in% c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)","C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)","C - 5 (COPEPODITE V)" ,
                                                                    "A + J (ADULT/JUVENILE)","JUVENILE","ADULT") ~ SPECIMEN_FORM %in% c("L", "H"),
 
## EUPUHASIIDS: T. RASCHII #juveniles -double check with ROMS on stages here 
    #Old 80's - 2010
    TAXON_NAME == "Thysanoessa raschii" & YEAR < 2011  ~ SPECIMEN_FORM %in% c("F","A") , 
    #Middle 2011-2018
    TAXON_NAME == "Thysanoessa raschii"  & YEAR > 2010 & YEAR < 2019 ~ SPECIMEN_FORM %in% c("F","A"), 
    #Recent 2019 - present
    TAXON_NAME == "Thysanoessa raschii"  & YEAR > 2018 ~ SPECIMEN_FORM %in% c("F", "K"), 
                                                                
## EUPUHASIIDS:  T. INERMIS
    #Old 80's - 2010
    TAXON_NAME == "Thysanoessa inermis" & YEAR < 2011  ~ SPECIMEN_FORM %in% c("F","A") , 
    #Middle 2011-2018
    TAXON_NAME == "Thysanoessa inermis"  & YEAR > 2010 & YEAR < 2019 ~ SPECIMEN_FORM %in% c("F","A"), 
    #Recent 2019 - present
    TAXON_NAME == "Thysanoessa inermis"  & YEAR > 2018 ~ SPECIMEN_FORM %in% c("F", "K"), 
    TRUE ~ FALSE))  %>%
    
  dplyr::select(-X) #remove extra column 
  

#### EVENTUALLY WILL INSERT WW DW AND C CONVERSIONS HERE!  #### 

#at some point we may want to change to ROMS names... later. 

  # group_by(ROMS_SP_NAME,CRUISE,YEAR,MONTH,DAY,GMT_DATE_TIME_TXT, HAUL_ID, HAUL_NAME, GEAR_NAME, MESH, 
  #          NET, HAUL_PERFORMANCE, LAT, LON, STATION_NAME, MAX_GEAR_DEPTH, MIN_GEAR_DEPTH, SAMPLE_DEPTH, BOTTOM_DEPTH, 
  #          DIS_PERVOLM2,DIS_PERVOLM3, SPECIMEN_FORM, STAGE_NAME) %>%
  # dplyr::summarise(EST_NUM_PERM2 = sum(EST_NUM_PERM2 ), EST_NUM_PERM3 = sum(EST_NUM_PERM3))

################################################################################################################
#make some plots to check the filtering before grouping across life stages 
################################################################################################################

df <- unique(zoop_filtered[c("YEAR", "TAXON_NAME", "GEAR_NAME", "MESH", "STAGE_NAME","SPECIMEN_FORM")])

unique(df$TAXON_NAME)

df<-df %>% 
  filter(!is.na(TAXON_NAME)) %>%
  filter(!STAGE_NAME == "NOT DETERMINED")
 
stage_plot<-df %>%
  mutate(STAGE_NAME = factor(STAGE_NAME, levels = c("C - 1 (COPEPODITE I)","C - 2 (COPEPODITE II)",
                                                       "C - 3 (COPEPODITE III)","C - 4 (COPEPODITE IV)",
                                                       "C - 5 (COPEPODITE V)","JUVENILE",  "ADULT", "A + J (ADULT/JUVENILE)"))) %>%
  ggplot( ) +
  geom_point(aes(x=YEAR, y = SPECIMEN_FORM, color = GEAR_NAME)) +
  facet_grid(TAXON_NAME~STAGE_NAME,labeller = label_wrap_gen(width=10))+
  theme_classic() +
  theme(axis.text.y = element_text(size=7),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.text.x = element_text(size = 5),
        strip.text.y = element_text(size = 5))
stage_plot

# pdf("output/plot_zoop_stages.pdf")
# stage_plot
# dev.off()

################################################################################################################
# Sum across lifestages  
################################################################################################################
zoop_grouped <- zoop_filtered %>%
                  group_by(TAXON_NAME, CRUISE,YEAR,MONTH,DAY,GMT_DATE_TIME_TXT, HAUL_ID, HAUL_NAME, GEAR_NAME, MESH, 
                            NET, HAUL_PERFORMANCE, LAT, LON, STATION_NAME, MAX_GEAR_DEPTH, MIN_GEAR_DEPTH, SAMPLE_DEPTH, BOTTOM_DEPTH, 
                            DIS_PERVOLM2,DIS_PERVOLM3) %>%
                  dplyr::summarise(EST_NUM_PERM2 = sum(EST_NUM_PERM2), EST_NUM_PERM3 = sum(EST_NUM_PERM3))

################################################################################################################
#expand the data set to add zero's 
################################################################################################################
#use this for all hauls to  add zeros into the data. 
Haul_Master <- unique(zoopdata[c("CRUISE", "DAY", "GEAR_NAME", "GMT_DATE_TIME_TXT", "LAT", "LON", "MAX_GEAR_DEPTH", 
                                 "MESH", "MONTH", "STATION_NAME",  "YEAR", "HAUL_ID", "BOTTOM_DEPTH")])

#list of all non-euphasiid species and stage names:
df_sp_list <- data.frame(ROMS_SP_NAME=sp_list)

#check to make sure all cruises are there
Haul_Master_Summary <- Haul_Master %>%
  group_by( CRUISE) %>%
  dplyr::summarise( n())

#create a df where all hauls exist for each species, left join on to this to add in all zeros
Master<-Haul_Master %>%
  expand_grid(sp_list) %>%
  dplyr::rename(TAXON_NAME = sp_list)

#join to get expanded hauls for each species hwich includes 0's
zoop_complete <- left_join(Master, zoop_grouped) %>%
  mutate(ROMS_SP_NAME  = case_when(TAXON_NAME %in% c("Calanus marshallae", "Calanus pacificus", "Calanus glacialis") ~ "Calanus_spp",
                                   TAXON_NAME %in% neocalanus_spp ~ "Neocalanus_spp",
                                   TAXON_NAME %in% pseudocalanus ~ "Pseudocalanus",
                                   TRUE  ~ TAXON_NAME)) %>%
  #change NA to 0's
  mutate(EST_NUM_PERM2=replace_na(EST_NUM_PERM2,0),EST_NUM_PERM3= replace_na(EST_NUM_PERM3,0)) %>%
  filter(!LAT >65, !LON > -154) # do rough spatial filter, will need to do this better later, mostly removing stuff in GOA

calanus <- zoop_complete %>%
  filter(TAXON_NAME %in% c("Calanus marshallae", "Calanus pacificus", "Calanus glacialis")) %>%
  filter(!YEAR < 2013)

#add on domain----

#Load shapefile that Dave sent for Ortiz regions
ortiz_shp <- st_read("data/ortiz_regions/BSIERP_regions_2012.shp")
ortiz_shp <- st_transform(ortiz_shp, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")#Tranform into WGS84 coordinate system
ortiz_shp <- subset(ortiz_shp, !DOMAIN == 15)
sf_ortiz = st_as_sf(ortiz_shp) #convert your shapefile to sf, too. 

sf_df = st_as_sf(calanus, coords = c("LON", "LAT"), crs = 4326) # change df to sf

#First do a join that contains the points to filter points that are tootally ouotside the regiono. 
#Then do a NN join to hopefully get a point or each polygon? 
df_contains <- st_join(sf_ortiz, sf_df, join = st_contains) 

#need to do some tidying to, rejoin OG lat longs and plot to look at trim
st_geometry(df_contains)<-NULL

df_filtered<-df_contains %>%
  filter(!name == "AK peninsula") %>%
  left_join(calanus) 
 
#write_csv(df_filtered, "data/Calanus_data_14-16.csv")
