library(sdmTMB)
library(dplyr)
###wd
setwd("~/Dropbox/GitHub/groundfish-o2-thresholds")

#Load functions
source("code/helper_funs.R")

###Oxygen data
##Combine with in situ data
insitu <- readRDS("data/processed_data/o2/insitu_combined.rds")
#Just columns of interest
dat <- insitu[,c("survey", "depth", "year", "date", "latitude", "longitude", "X", "Y", "temperature_C", "do_mlpL", "salinity_psu", "sigma0_kgm3", "O2_umolkg", "event_id", "doy")]
#Drop if NA in O2
dat <- drop_na(dat, latitude, longitude, O2_umolkg, temperature_C, salinity_psu, year)
#Calculate pO2
#Calculate pO2 from umol kg
dat$po2 <- calc_po2_sat(salinity=dat$salinity_psu, temp=dat$temperature_C, depth=dat$depth, oxygen=dat$O2_umolkg, lat=dat$latitude, long=dat$longitude, umol_m3=T, ml_L=F)

##Add region
# load regional polygons
regions.hull <- readRDS("data/processed_data/regions_hull.rds")
#make dataframe an sf object
dat_df <-  st_as_sf(dat, coords = c("longitude", "latitude"), crs = st_crs(4326))
dat_df$latitude <- dat$latitude
dat_df$longitude <- dat$longitude
# cycle through all regions
region_list <- c("ai", "bc", "cc", "ebs", "goa")
dats <- list()
for (i in 1:length(region_list)) {
  region <- region_list[i]
  poly <- regions.hull[i,2]
  # pull out observations within each region
  region_dat  <- st_filter(dat_df, poly)
  region_dat$region <- paste(region_list[i])
  dats[[i]] <- as.data.frame(region_dat)
}
#Bind back together
dat <- bind_rows(dats)
dat <- unique(dat)

#Remove weird depths
dat <- filter(dat, depth>0)

#Remove oxygen outliers
dat <- filter(dat, O2_umolkg<1500)

#Remove ai from list
region_list <- region_list[-1]

#Log depth
dat$depth_ln <- log(dat$depth)

##Get grids set up
##Load grids
load("data/nwfsc_grid.rda")
load("data/dfo_synoptic_grid.rda")

##Make a list of grids
grids <- list(dfo_synoptic_grid, nwfsc_grid)

#Fit models for sigma0, temperature, and O2 and predict to grid
preds <- list()
for(i in 1:length(region_list)){
  dat.2.use <- filter(dat, region==region_list[i])
  dat.2.use$year <- as.numeric(dat.2.use$year)
  extra_years <- c(seq(2008,2024,1))
  #Create mesh
  ## Make Mesh and fit model ####
  spde <- make_mesh(data = dat.2.use,
                    xy_cols = c("X", "Y"),
                    cutoff = 45)
  #Fit model
  formula <- po2 ~ 1+s(depth_ln)+s(doy)
  
    m1 <- try(sdmTMB(formula = as.formula(formula),
    mesh = spde,
    data = dat.2.use,
    family = gaussian(),
    time = "year",
    spatial = "on",
    spatiotemporal  = "ar1",
    extra_time=c(extra_years)
    ))

    formula <- temperature_C ~ 1+s(depth_ln)+s(doy)
    
    m2 <- try(sdmTMB(formula = as.formula(formula),
                     mesh = spde,
                     data = dat.2.use,
                     family = gaussian(),
                     time = "year",
                     spatial = "on",
                     spatiotemporal  = "ar1",
                     extra_time=c(extra_years)
    ))
    
    formula <- sigma0_kgm3 ~ 1+s(depth_ln)+s(doy)
    m3 <- try(sdmTMB(formula = as.formula(formula),
                     mesh = spde,
                     data = dat.2.use,
                     family = gaussian(),
                     time = "year",
                     spatial = "on",
                     spatiotemporal  = "ar1",
                     extra_time=c(extra_years)
    ))

    grid <- grids[[i]]
    region.2.use <- region_list[i]
    
    #Remove negative depths
    grid <- filter(grid, depth_m>0)
    #Make depth ln
    grid$depth_ln <- log(grid$depth_m)
    #Get X and Y
    grid$longitude <- grid$lon
    grid$latitude <- grid$lat
    # UTM transformation
    dat_ll = grid
    sp::coordinates(dat_ll) <- c("lon", "lat")
    sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
    # convert to utm with spTransform
    dat_utm = sp::spTransform(dat_ll, 
                              sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
    # convert back from sp object to data frame
    grid = as.data.frame(dat_utm)
    grid = dplyr::rename(grid, X = coords.x1,
                         Y = coords.x2)
    #Expand by years
    grid <- grid %>%
      expand_grid(year = extra_years)
    grid$year <- as.numeric(grid$year)
    
    #Add doy
    #July 1st
    grid$doy <- 182
    
    #Add region column
    grid$region <- region_list[i]
    #Predict pO2
    p <- predict(m1, newdata = grid)
    grid <- grids[[i]]
  m <- region_models[[i]]
  region.2.use <- region_list[i]
  dat.2.use <- filter(dat, region==region.2.use)
  print(region_list[i])
  #Remove negative depths
  grid <- filter(grid, depth_m>0)
  #Make depth ln
  grid$depth_ln <- log(grid$depth_m)
  #Get X and Y
  grid$longitude <- grid$lon
  grid$latitute <- grid$lat
  # UTM transformation
  dat_ll = grid
  sp::coordinates(dat_ll) <- c("lon", "lat")
  sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = sp::spTransform(dat_ll, 
                            sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  grid = as.data.frame(dat_utm)
  grid = dplyr::rename(grid, X = coords.x1,
                       Y = coords.x2)
  #Expand by years
  grid <- grid %>%
    expand_grid(year = extra_years)
  grid$year <- as.numeric(grid$year)
  
  #Add doy
  #July 1st
  grid$doy <- 182
  
  #Add region column
  grid$region <- region_list[i]
  #Predict pO2
  p1 <- predict(m1, newdata = grid)
  p2 <- predict(m2, newdata = grid)
  p3 <- predict(m3, newdata = grid)

  grid$po2 <- p1$est
  grid$temperature_C <- p2$est
  grid$sigma0_kgm3 <- p3$est
  
  preds[[i]] <- grid
}


saveRDS(preds[[1]], "data/processed_data/o2/bc_predictions_grid.rds")
saveRDS(preds[[2]], "data/processed_data/o2/cc_predictions_grid.rds")
