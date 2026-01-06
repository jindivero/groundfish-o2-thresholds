install.packages("pak")
library(pak)
pak::pkg_install("DFO-NOAA-Pacific/surveyjoin")
library(surveyjoin)
library(dplyr)
library(tidyr)
library(readxl)
library(sf)

#Set wd
setwd("~/Dropbox/GitHub/groundfish-o2-thresholds")

##Initial data cache
cache_data()
load_sql_data()

#All data
dat <- get_data()

#Load functions
source("code/helper_funs.R")

#Rename troublesome common names
dat$common_name <- ifelse(dat$common_name=="pacific spiny dogfish", "spiny dogfish", dat$common_name)
dat$common_name <- ifelse(grepl("rougheye", dat$common_name), "rougheye rockfish", dat$common_name)
#Rename co
dat$latitude <- dat$lat_start
dat$longitude <- dat$lon_start
dat$depth <- dat$depth_m
dat$depth_m <- NULL

dat <- dat  %>%
  drop_na(latitude, longitude)

# UTM transformation
dat_ll = dat
sp::coordinates(dat_ll) <- c("longitude", "latitude")
sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
# convert to utm with spTransform
dat_utm = sp::spTransform(dat_ll, 
                          sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
# convert back from sp object to data frame
dat = as.data.frame(dat_utm)
dat = dplyr::rename(dat, X = coords.x1,
                    Y = coords.x2)
dat$latitude <- dat$lat_start
dat$longitude <- dat$lon_start

#Make date POSIXct
dat$date <- as.POSIXct(dat$date)

#Some other columns
dat$survey <- case_when(grepl("SYN", dat$survey_name)~"dfo",
                        (grepl("NWFSC", dat$survey_name)|dat$survey_name=="Triennial")~"nwfsc",
                        dat$survey_name=="Aleutian Islands"|dat$survey_name=="eastern Bering Sea"|dat$survey_name=="northern Bering Sea"|dat$survey_name=="Bering Sea Slope"~"afsc_bsai",
                        dat$survey_name=="Gulf of Alaska"~"afsc_goa")

dat$region<- case_when(dat$survey_name=="northern Bering Sea"~"nbs",
                       dat$survey_name=="eastern Bering Sea"|dat$survey_name=="Bering Sea Slope"~"ebs",
                       dat$survey=="nwfsc"~"cc",
                       dat$survey=="afsc_goa"~"goa",
                       dat$survey=="dfo"~"bc")

#Combine with in situ data
insitu <- readRDS("data/processed_data/o2/insitu_combined.rds")
#Remove IPHC, because this is already in the IPHC fish catch data
insitu <- filter(insitu, survey!="iphc")
insitu <- unique(insitu)

#Just columns of interest
dat <- left_join(dat, insitu[,c("temperature_C", "do_mlpL", "salinity_psu", "sigma0_kgm3", "O2_umolkg", "event_id")], by="event_id")

#Note: current taxa options for calculating metabolic index
species <- read_excel("data/species_table.xlsx")
species$common_name <- tolower(species$common_name)
species$scientific_name <- tolower(species$scientific_name)

#Ones w/ IPHC data
sub_species <- c("sablefish", "pacific cod", "pacific halibut", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")
sub <-filter(species, common_name %in% sub_species)

spcs <- tolower(sub$common_name)
sci_names <- tolower(sub$scientific_name)
taxas <- tolower(sub$MI_Taxa)
file_names <- spcs

#Do the full MI equation (fancy=T) or the abbreviated (just Eo; fancy=F)
fancy <- F
#Ones to include IPHC data
for(i in 1:length(spcs)){
  spc <- spcs[i]
  print(spc)
  sci_name <- sci_names[i]
  file_name <- file_names[i]
  taxa <- taxas[i]
  try(prepare_data2(dat, spc=spc, sci_name=sci_name,  taxa=taxa, iphc=T, file_name=file_name,fancy=fancy))
  gc()
}

#To not include IPHC data
sub_species <- c("sablefish", "pacific cod", "pacific halibut", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")
sub <-filter(species, !(common_name %in% sub_species))
spcs <- tolower(sub$common_name)
sci_names <- tolower(sub$scientific_name)
taxas <- tolower(sub$MI_Taxa)
file_names <- spcs
#Ones to include IPHC data
for(i in 1:length(spcs)){
  spc <- spcs[i]
  print(spc)
  sci_name <- sci_names[i]
  file_name <- file_names[i]
  taxa <- taxas[i]
  try(prepare_data2(full_data=dat, spc=spc, sci_name=sci_name, taxa=taxa, iphc=F, file_name=file_name, fancy))
  gc()
  
}

