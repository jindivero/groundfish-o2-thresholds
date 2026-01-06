remove.packages("sdmTMB")
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,  ref="newbreakpt")
library(sdmTMB)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)

###wd
setwd("~/Dropbox/GitHub/groundfish-o2-thresholds")

#Load functions
source("code/helper_funs.R")

#Output folder
output_folder <- "region_comp4"

#ggplot themes
theme_set(theme_bw(base_size = 18))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

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

#Add coastwide dataset
coastwide <- dat
coastwide$region <- "coastwide"
dat <- bind_rows(dat,coastwide)

#Plot O2 by depth and temp in each region
dat <- filter(dat, region!="ai")
dat$region <- factor(dat$region, levels=c("ebs", "goa", "bc", "cc", "coastwide"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current", "Coastwide")
names(labs) <- c("ebs", "goa", "bc", "cc", "coastwide")
ggplot(filter(dat, region!="coastwide"), aes(x=po2, y=depth))+
  geom_point(aes(colour=temperature_C), size=0.5, alpha=0.1)+
  facet_wrap("region", labeller = as_labeller(labs), scales="free_x")+
  scale_colour_viridis(option="inferno", name="Temperature (C)")+
  scale_y_reverse()+
  xlab("Partial Pressure Oxygen (kPa)")+
  ylab("Depth (m)")+
  theme(legend.position = "top")

ggsave(
  paste0("output/", output_folder, "/plots_final/O2_depth_temp.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height= 8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#Correlation matrix
library(purrr)

region_list <- c("cc", "bc", "goa", "ebs", "coastwide")

for(i in 1:length(region_list)){
  dat.2.use <- filter(dat, region==region_list[i])
  cor_mat <- cor(select(dat.2.use, temperature_C, depth, po2), use = "pairwise.complete.obs")
  cor_long <- as.data.frame(as.table(cor_mat))
  cor_long$region <- region_list[i]
  if(i==1){
    cor_all <- cor_long
  } else {
    cor_all <- bind_rows(cor_all, cor_long)
  }
  #Remove duplicates
  cor_all <- cor_all %>%
    filter(as.character(Var1) > as.character(Var2))
}

#Plot correlogram
#Labels
cor_all$region <- factor(cor_all$region, levels=c("ebs", "goa", "bc", "cc", "coastwide"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current", "Coastwide")
names(labs) <- c("ebs", "goa", "bc", "cc", "coastwide")
#replace temperature_C with Temp
cor_all$Var1 <- gsub("temperature_C", "Temp", cor_all$Var1)
cor_all$Var2 <- gsub("temperature_C", "Temp", cor_all$Var2)
cor_all$Var1 <- gsub("depth", "Depth", cor_all$Var1)
cor_all$Var2 <- gsub("depth", "Depth", cor_all$Var2)
cor_all$Var1 <- gsub("po2", "pO2", cor_all$Var1)
cor_all$Var2 <- gsub("po2", "pO2", cor_all$Var2)

ggplot(cor_all, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = round(Freq, 2)), size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, limits = c(-1, 1)) +
  facet_wrap(~region, labeller=labeller(region=labs)) +
  coord_fixed() +
  labs(x = NULL, y = NULL, fill = "Correlation")+
  theme_minimal(base_size=17)+
  theme(legend.position=c(0.8,0.2))

#Save
ggsave(
  paste0("output/", output_folder, "/plots/O2_depth_temp_correlelogram.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height= 8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)
#Fit model just to year (each region sep)
year_models <- list()
for(i in 1:length(region_list)){
  dat.2.use <- filter(dat, region==region_list[i])
  dat.2.use$year <- as.factor(dat.2.use$year)
  #Create mesh
  ## Make Mesh and fit model ####
  spde <- make_mesh(data = dat.2.use,
                    xy_cols = c("X", "Y"),
                    cutoff = 45)
  #Fit model
    formula <- po2 ~ 1+as.factor(year)+s(depth_ln)+s(doy)
  
  m1 <- try(sdmTMB(
    formula = as.formula(formula),
    mesh = spde,
    time="year",
    data = dat.2.use,
    family = gaussian(),
    spatial = "on",
    spatiotemporal  = "iid"
  ))
  year_models[[i]] <- m1
}

#What are the annual effects?
for(i in 1:length(year_models)){
  m <- year_models[[i]]
  #Get the coefficients
  summary(m)
  #Create a dataframe of the coefficients
  if(i==1){
    coefs <- data.frame(coef(m))
    colnames(coefs) <- region_list[i]
    coefs$term <- row.names(coefs)
    coefs <- select(coefs, term, bc)
  }
  if(i>1){
    coefs.temp <- data.frame(coef(m))
    colnames(coefs.temp) <- region_list[i]
    coefs.temp$term <- row.names(coefs.temp)
    coefs <- left_join(coefs, coefs.temp, by="term")
  }
}

#What is the SD of spatial & spatiotmeporal variation
for(i in 1:length(year_models)){
  m <- year_models[[i]]
  #Get the coefficients
  a <- tidy(m, effects = "ran_pars")
  a <- a[,1:2]
  colnames(a) <- c("term", paste(region_list[i]))
  #Add to the year terms
  if(i==1){
  pars <- a
  } else {
  pars <- left_join(pars, a, by="term")
  }
}


###Predict to grid
##Load grids
load("data/nwfsc_grid.rda")
load("data/afsc_grid.rda")
ebs_grid <- filter(afsc_grid, survey=="Eastern Bering Sea Crab/Groundfish Bottom Trawl Survey"|survey=="Eastern Bering Sea Slope Bottom Trawl Survey")       
goa_grid <- filter(afsc_grid, survey=="Gulf of Alaska Bottom Trawl Survey")
load("data/dfo_synoptic_grid.rda")

###Add depth to ebs and goa grids from NOAA bathymetry data
##Bathymetry data
# load bathymetry data
bathy_all <- readRDS("data/bathymetry_regions_from_grids.rds")
#Just goa and ebs
bathy_all <- bathy_all[4:5]

# fit models for each region
make_depth_model <- function(bathydat) {
  
  spde <- make_mesh(data = as.data.frame(bathydat), xy_cols = c("X", "Y"), n_knots = 300)
  
  depth_model <- sdmTMB(log(noaadepth) ~ 1,
                        data = as.data.frame(bathydat),
                        spatial = "on", 
                        mesh = spde,
                        family = gaussian()
  )
}
depth_models <- lapply(X = bathy_all,
                       FUN = make_depth_model)


# load all data
# cycle through the alaska grids
afsc_grids <- list(ebs_grid, goa_grid)

for (i in 1:length(afsc_grids)) {
  grid2use <- afsc_grids[[i]]
  #Get X and Y
  grid2use$longitude <- grid2use$lon
  grid2use$latitude <- grid2use$lat
  # UTM transformation
  dat_ll = grid2use
  sp::coordinates(dat_ll) <- c("lon", "lat")
  sp::proj4string(dat_ll) <- sp::CRS("+proj=longlat +datum=WGS84")
  # convert to utm with spTransform
  dat_utm = sp::spTransform(dat_ll, 
                            sp::CRS("+proj=utm +zone=10 +datum=WGS84 +units=km"))
  # convert back from sp object to data frame
  dat3 = as.data.frame(dat_utm)
  dat3= dplyr::rename(dat3, X = coords.x1,
                       Y = coords.x2)
  dat3$lat <- dat3$latitude
  dat3$lon <- dat3$longitude
  model.2.use <- depth_models[[i]]
  
  # get predicted log(depth) for each observation, based on model fit to that region
  region_dat_predict <- predict(model.2.use, as.data.frame(dat3))
  # add the predicted log(depth) to the data frame
  region_dat_predict$depth_m <- exp(region_dat_predict$est)
  dat2 <- select(region_dat_predict, area, depth_m, survey, survey_domain_year, lon, lat)
  
  if(i==1){
  ebs_grid <- dat2
  } else {
  goa_grid <- dat2
  }
}

##Make a list of grids
grids <- list(dfo_synoptic_grid, nwfsc_grid, ebs_grid, goa_grid)

##Predict pO2 across grid
preds <- list()
for(i in 1:length(region_list)){
  grid <- grids[[i]]
  m <- year_models[[i]]
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
    expand_grid(year = unique(dat.2.use$year)) %>%
    mutate(year = as.factor(year))
  grid$year <- as.factor(grid$year)
  
  #Add doy
  #July 1st
  grid$doy <- 182
  
  #Add region column
  grid$region <- region_list[i]
  #Predict pO2
  p <- predict(m, newdata = grid)
  #Add to list
  preds[[i]] <- p
}

#Bind together
preds <- bind_rows(preds)

#Map spatial variation in each region
###Map#########
map_data <- rnaturalearth::ne_countries(scale = "large",
                                        returnclass = "sf",
                                        continent = "North America")

us_coast_proj <- sf::st_transform(map_data, crs = 32610)

for(i in 1:length(region_list)){
  grid2use<- filter(preds, region==region_list[i])
  grid2use<- filter(grid2use, year=="2021")
  print(region_list[i])

ggplot(us_coast_proj) + geom_sf() +
    geom_point(grid2use, mapping=aes(x=X*1000, y=Y*1000,colour=omega_s), size=0.1)+
    xlim(min(grid2use$X)*1000, max(grid2use$X)*1000)+
    ylim(min(grid2use$Y)*1000, max(grid2use$Y)*1000)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(axis.text.x=element_blank())+
    scale_colour_viridis()

ggsave(
  paste("output/plots/o2_omega_",region_list[i], ".png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

}

##Spatiotemporal effects for each year
for(i in 1:length(region_list)){
  grid2use<- filter(preds, region==region_list[i])
  print(region_list[i])
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(grid2use, mapping=aes(x=X*1000, y=Y*1000,colour=epsilon_st), size=0.1)+
    xlim(min(grid2use$X)*1000, max(grid2use$X)*1000)+
    ylim(min(grid2use$Y)*1000, max(grid2use$Y)*1000)+
      facet_wrap("year", ncol=4)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(axis.text.x=element_blank())+
    scale_colour_viridis()
  
  ggsave(
    paste("output/plots/o2_epsilon_",region_list[i], ".png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height = 11,
    units = c("in"),
    dpi = 600,
    limitsize = TRUE, bg="white"
  )
  
}

#Calculate the coefficient of variation of epsilon_st for each latitude & longitude combination
epsilon_cv <- preds %>%
  group_by(X, Y, region) %>%
  summarise(epsilon_st_sd = sd(epsilon_st),epsilon_st_mean = mean(epsilon_st),epsilon_st_var = var(epsilon_st), epsilon_st_cv = (sd(epsilon_st)/mean(epsilon_st))) %>%
  ungroup()

#plot
for(i in 1:length(region_list)){
  grid2use<- filter(epsilon_cv, region==region_list[i])
  print(region_list[i])
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(grid2use, mapping=aes(x=X*1000, y=Y*1000,colour=epsilon_st_var), size=0.1)+
    xlim(min(grid2use$X)*1000, max(grid2use$X)*1000)+
    ylim(min(grid2use$Y)*1000, max(grid2use$Y)*1000)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(axis.text.x=element_blank())+
    scale_colour_viridis()
  
  ggsave(
    paste("output/plots/o2_epsilon_cv_",region_list[i], ".png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height = 5,
    units = c("in"),
    dpi = 600,
    limitsize = TRUE, bg="white"
  )
  
}

##Proportion of annual, spatial, and spatio-temporal variation
##Add value of column in coefs matching region and year of preds
#Add year column
#last four digits of coefs$term
coefs$year <- substr(coefs$term, nchar(coefs$term)-3, nchar(coefs$term))
#pivot long
coefs2 <- pivot_longer(coefs, cols=2:5, names_to="region", values_to="value")
#Match to preds
preds <- left_join(preds, coefs2, by=c("region", "year"))

#Rename
preds$annual_fixed <- preds$value

#Calculate proportion of prediction from each term
preds$spatial <- abs(preds$omega_s)/(abs(preds$omega_s)+abs(preds$epsilon_st)+abs(preds$annual_fixed))
preds$spatiotemp <- abs(preds$epsilon_st)/(abs(preds$omega_s)+abs(preds$epsilon_st)+abs(preds$annual_fixed))
preds$annual <- abs(preds$annual_fixed)/(abs(preds$omega_s)+abs(preds$epsilon_st)+abs(preds$annual_fixed))

###These look terrible
#Ternary plot
install.packages("ggtern")
library(ggtern)

ggtern(data=filter(preds,year=="2021"),aes(annual,spatial,spatiotemp)) + 
geom_mask() +
  facet_wrap("region")+
  geom_point(size=0.2, alpha=0.5, aes(colour=region)) + 
  theme_bw() +
  theme_showarrows() +
  theme_clockwise()+
  ggtitle("2021")

ggtern(data=filter(preds,year=="2012"),aes(annual,spatial,spatiotemp)) + 
  geom_mask() +
  facet_wrap("region")+
  geom_point(size=0.2, alpha=0.5, aes(colour=region)) + 
  theme_bw() +
  theme_showarrows() +
  theme_clockwise()+
  ggtitle("2012")

library(lattice)
cloud(abs(epsilon_st)~ abs(omega_s) * abs(annual_fixed), pch = ".", data = filter(preds, region=="cc"), group=year)
cloud(abs(epsilon_st)~ abs(omega_s) * abs(annual_fixed), pch = ".", data = filter(preds, region=="bc"), group=year)
cloud(abs(epsilon_st)~ abs(omega_s) * abs(annual_fixed), pch = ".", data = filter(preds, region=="goa"), group=year)
cloud(abs(epsilon_st)~ abs(omega_s) * abs(annual_fixed), pch = ".", data = filter(preds, region=="ebs"), group=year)

#Overall
pars2 <- bind_rows(coefs, pars)
pars2 <- filter(pars2, term!="(Intercept)"&term!="range"&term!="phi")
#Pivot long
pars2 <- pivot_longer(pars2, cols=2:5, names_to="region", values_to="value")
pars2$type <- ifelse(grepl("year", pars2$term), "fixed", "random")
#Plot
ggplot(pars2, aes(x=region, y=abs(value)))+
  geom_col(aes(fill=term, alpha=type))+
  theme_minimal(base_size=12)+
  scale_alpha_manual(values=c(0.5, 1))

#Plot spatial marginal effects
#plot
preds$spatial_marg_prop <- (preds$est-preds$omega_s)/preds$est
preds$spatiotemp_marg_prop <- (preds$est-preds$epsilon_st)/preds$est

for(i in 1:length(region_list)){
  grid2use<- filter(preds, region==region_list[i])
  print(region_list[i])
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(filter(grid2use,est>0), mapping=aes(x=X*1000, y=Y*1000,colour=spatial_marg_prop), size=0.1)+
    xlim(min(grid2use$X)*1000, max(grid2use$X)*1000)+
    ylim(min(grid2use$Y)*1000, max(grid2use$Y)*1000)+
    theme_minimal(base_size=12)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(axis.text.x=element_blank())+
    scale_colour_viridis()
  
  ggsave(
    paste("output/plots/o2_spatial_prop_change",region_list[i], ".png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height = 5,
    units = c("in"),
    dpi = 600,
    limitsize = TRUE, bg="white"
  )
  
}
