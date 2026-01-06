library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(openxlsx2)
library(ggpubr)
install.packages("concaveman")
library(concaveman)
library(sf)

setwd("~/Dropbox/GitHub/groundfish-o2-thresholds")

#Load functions
source("code/helper_funs.R")

#Output folder
output_folder <- "region_comp4"

#Plot themes
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

##Set up mapping
map_data <- rnaturalearth::ne_countries(scale = "large",
                                        returnclass = "sf",
                                        continent = "North America")

us_coast_proj <- sf::st_transform(map_data, crs = 32610)

#Load data
dat <- list.files(path = "data/processed_data/fish3", pattern = ".rds", full.names=T) %>%
  map(readRDS) %>% 
  bind_rows()

#Remove Aleutian Islands and Northern Bering Sea
dat <- filter(dat, region!="ai")
dat <- filter(dat, region!="nbs")

#Combine IPHC and bottom trawl catch data
dat$catch_weight_combined <- ifelse(is.na(dat$catch_weight), dat$cpue_weight, dat$catch_weight)
dat$catch_count_combined <- ifelse(is.na(dat$catch_numbers), dat$cpue_count, dat$catch_numbers)

#Clean up
dat$catch_weight_combined <- replace(dat$catch_weight_combined, dat$catch_weight_combined == "Inf", NA)
dat$catch_count_combined <- replace(dat$catch_count_combined, dat$catch_count_combined == "Inf", NA)
dat$catch_weight_combined <- replace(dat$catch_weight_combined, dat$catch_weight_combined == "NaN", NA)
dat$catch_count_combined <- replace(dat$catch_count_combined, dat$catch_count_combined == "NaN", NA)

##Get depth habitat for filtering
species_table <- read_excel("data/species_table.xlsx")
species_table$common_name <- tolower(species_table$common_name)
species <- unique(species_table$common_name)
#IPHC species that use counts and not weights
species_iphc <- c("sablefish", "pacific cod", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")

bottom_trawl_only <- F

for(i in 1:length(species)){
  this_species <- species[i]
  dat_depth2 <- filter(dat, common_name==this_species)
  print(this_species)
  
  # Sort by depth
  dat_depth2 <- dat_depth2[order(dat_depth2$depth), ]
  
  #Calculate the cumulative sum of catch by depth
  if(this_species %in% species_iphc){
    dat_depth2$catch2use <- dat_depth2$catch_count_combined
  } else {
    dat_depth2$catch2use <- dat_depth2$catch_weight_combined
  } 
  if(bottom_trawl_only){
    dat_depth2$catch2use <- dat_depth2$catch_weight
  }
  
  dat_depth2 <- drop_na(dat_depth2, catch2use)
  dat_depth2$cumsum_catch <- cumsum(dat_depth2$catch2use)

  #Calculate the proportional cumulative sum
  dat_depth2$prop_cumsum_catch <- dat_depth2$cumsum_catch / sum(dat_depth2$catch2use, na.rm=T)
  
  # Find the index of the closest value to 99% (0.99) in prop_cumsum_var1
  closest_index <- which.min(abs(dat_depth2$prop_cumsum_catch - 0.99))
  
  # Get the depth at 99% catch
  closest_value <- dat_depth2$depth[closest_index]
  
  #Add to dataframe
  dat_depth2$filtered_depth <- closest_value
  
  #Dataframe
  depth2use <- data.frame(common_name=this_species, depth=closest_value)
  
  if(i==1){
    depths <- depth2use
    dat_depths <- dat_depth2
  } else {
    depths <- bind_rows(depths, depth2use)
    dat_depths <- bind_rows(dat_depths, dat_depth2)
  }
}

#Add depths to species table
species_table <- left_join(species_table, depths, by="common_name")
write_xlsx(species_table, "data/species_table.xlsx")

#Plot cumulative sum by depth
ggplot(dat_depths, aes(x=depth, y=prop_cumsum_catch))+
  geom_line()+
  geom_vline(species_table, mapping=aes(xintercept=depth), linetype="dashed")+
  facet_wrap("common_name", ncol=4, labeller=labeller(common_name=label_wrap_gen(15)))+
  ylab("Cumulative Sum of Catch")+
  xlab("Depth (m)")

ggsave(
  paste0("output/", output_folder, "/plots/cumulative_depth.png"),
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

###Filter geographically--latitude and longitude
###Filter for latitude
#Latitude
for(i in 1:length(species)){
  this_species <- species[i]
  dat2 <- filter(dat, common_name==this_species)
  print(this_species)
  
  #Select catch data to use
  if(this_species %in% species_iphc){
    dat2$catch2use <- dat2$catch_count_combined
  } else {
    dat2$catch2use <- dat2$catch_weight_combined
  } 
  if(bottom_trawl_only){
    dat2$catch2use <- dat2$catch_weight
  }
  
  #Calculate the cumulative sum of catch by latitude in each region
  #Southern extent:
  dat_cc<- filter(dat2, region=="cc")
  dat_cc2 <- filter(dat_cc, catch2use>0)
  if(nrow(dat_cc2>50)){
  
  # Sort by latitude
  dat_cc <- dat_cc[order(dat_cc$latitude, decreasing=T), ]
  
  #drop NA
  dat_cc <- drop_na(dat_cc, catch2use)
  
  #Calculate cumulative sum
  dat_cc$cumsum_catch <- cumsum(dat_cc$catch2use)
  
  #Calculate the proportional cumulative sum
  dat_cc$prop_cumsum_catch <- dat_cc$cumsum_catch / sum(dat_cc$catch2use, na.rm=T)
  
  # Find the index of the closest value to 99% (0.99) in prop_cumsum_var1
  closest_index <- which.min(abs(dat_cc$prop_cumsum_catch - 0.99))
  
  # Get the depth at 99% catch
  closest_value <- dat_cc$latitude[closest_index]
  
  #Add to dataframe
  dat_cc$type <- "southern_limit"
  if(i==1){
    dat_limits <- dat_cc
    species_table$southern_limit[i] <- closest_value
  } else {
    dat_limits <- bind_rows(dat_limits, dat_cc)
    species_table$southern_limit[i] <- closest_value
  }
  }
  if(nrow(dat_cc2)<=50){
    dat_bc<- filter(dat2, region=="bc")
    dat_bc2 <- filter(dat_bc, catch2use>0)
    if(nrow(dat_bc2)>50){
      # Sort by latitude
      dat_bc <- dat_bc[order(dat_bc$latitude, decreasing=T), ]
      
      #drop NA
      dat_bc <- drop_na(dat_bc, catch2use)
      
      #Calculate cumulative sum
      dat_bc$cumsum_catch <- cumsum(dat_bc$catch2use)
      
      #Calculate the proportional cumulative sum
      dat_bc$prop_cumsum_catch <- dat_bc$cumsum_catch / sum(dat_bc$catch2use, na.rm=T)
      
      # Find the index of the closest value to 99% (0.99) in prop_cumsum_var1
      closest_index <- which.min(abs(dat_bc$prop_cumsum_catch - 0.99))
      
      # Get the depth at 99% catch
      closest_value <- dat_bc$latitude[closest_index]
      
      #Add to dataframe
      dat_bc$type <- "southern_limit"
      if(i==1){
        dat_limits <- dat_bc
        species_table$southern_limit[i] <- closest_value
      } else {
        dat_limits <- bind_rows(dat_limits, dat_bc)
        species_table$southern_limit[i] <- closest_value
      }
    }
  }
  
  ##Northern extent
  #If present in EBS
  dat_ebs <- filter(dat2, region=="ebs")
  dat_ebs2 <- filter(dat_ebs, catch2use>0)
  if(nrow(dat_ebs2>50)){
    # Sort by latitude
    dat_ebs <- dat_ebs[order(dat_ebs$latitude, decreasing=F),]
    #drop NA
    dat_ebs <- drop_na(dat_ebs, catch2use)
    #Calculate cumulative sum
    dat_ebs$cumsum_catch <- cumsum(dat_ebs$catch2use)
    
    #Calculate the proportional cumulative sum
    dat_ebs$prop_cumsum_catch <- dat_ebs$cumsum_catch / sum(dat_ebs$catch2use, na.rm=T)
    
    # Find the index of the closest value to 99% (0.99) in prop_cumsum_var1
    closest_index <- which.min(abs(dat_ebs$prop_cumsum_catch - 0.99))
    
    # Get the depth at 99% catch
    closest_value <- dat_ebs$latitude[closest_index]
    
    #Add to dataframe for i row
    species_table$ebs_limit[i] <- closest_value
    
    #Label dataset
    dat_ebs$type <- "ebs_limit"
    dat_limits <- bind_rows(dat_limits, dat_ebs)
  }
  if(nrow(dat_ebs2)<=50){
    #If not in EBS, check GOA
    species_table$ebs_limit[i] <- NA
    dat_goa <- filter(dat2, region=="goa")
    dat_goa2 <- filter(dat_goa, catch2use>0)
    if(nrow(dat_goa2>50)){
    # Sort by longitude (because trimming eastern extent, not by latitude)
    dat_goa <- dat_goa[order(dat_goa$longitude, decreasing=T), ]

    #drop NA
    dat_goa <- drop_na(dat_goa, catch2use)
    
    #Calculate cumulative sum
    dat_goa$cumsum_catch <- cumsum(dat_goa$catch2use)
    
    #Calculate the proportional cumulative sum
    dat_goa$prop_cumsum_catch <- dat_goa$cumsum_catch / sum(dat_goa$catch2use, na.rm=T)
    
    # Find the index of the closest value to 99% (0.99) in prop_cumsum_var1
    closest_index <- which.min(abs(dat_goa$prop_cumsum_catch - 0.99))
    
    # Get the depth at 99% catch
    closest_value <- dat_goa$longitude[closest_index]
    
    #Label dataframe
    dat_goa$type <- "goa_limit"
    
    #Add to dataframe for i row
    species_table$goa_limit[i] <- closest_value
    dat_limits <- bind_rows(dat_limits, dat_goa)
}
  }
    #If not enough in GOA, check in BC
    if(nrow(dat_goa2)<=50){
      dat_bc <- filter(dat2, region=="bc")
      dat_bc2 <- filter(dat_bc, catch2use>0)
      if(nrow(dat_bc2>50)){
        # Sort by latitude
        dat_bc <- dat_bc[order(dat_bc$latitude, decreasing=T), ]
        
        #drop NA
        dat_bc <- drop_na(dat_bc, catch2use)
        
        #Calculate cumulative sum
        dat_bc$cumsum_catch <- cumsum(dat_bc$catch2use)
        
        #Calculate the proportional cumulative sum
        dat_bc$prop_cumsum_catch <- dat_bc$cumsum_catch / sum(dat_bc$catch2use, na.rm=T)
        
        # Find the index of the closest value to 99% (0.99) in prop_cumsum_var1
        closest_index <- which.min(abs(dat_bc$prop_cumsum_catch - 0.99))
        
        # Get the depth at 99% catch
        closest_value <- dat_bc$latitude[closest_index]
        
        #Label dataset
        dat_bc$type <- "bc_limit"
        
        #Add to dataframe
        species_table$bc_limit[i] <- closest_value
        dat_lats2 <- bind_rows(dat_lats2, dat_bc)
    }
    }
}

#Check constraints--see how many positive catches are removed, plot, make sure reasonable
for(i in 1:nrow(species_table)){
  #Filter data for species
  this_species <- species_table$common_name[i]
  dat2 <- filter(dat, common_name==this_species)
  dat_pos <- filter(dat2, catch_weight_combined>0 | catch_count_combined>0)
  #Number of positive catches in complete dataset
  n_total <- nrow(filter(dat2, catch_weight_combined>0 | catch_count_combined>0))
  
  #Filter data to geographic coordinates
  dat_filtered <- filter(dat2, latitude >= species_table$southern_limit[i])
  if(!is.na(species_table$ebs_limit[i])){
    #Separate out the Bering Sea only--because Gulf of Alaska messes this up
    dat_filtered_ebs <- filter(dat_filtered, region=="ebs")
    dat_filtered <- filter(dat_filtered, region!="ebs")
    dat_filtered_ebs <- filter(dat_filtered_ebs, latitude <= species_table$ebs_limit[i])
    dat_filtered <- bind_rows(dat_filtered, dat_filtered_ebs)
  }
  if(!is.na(species_table$goa_limit[i])){
    dat_filtered <- filter(dat_filtered, longitude >= species_table$goa_limit[i])
  }
  
  #Number of positive catches in filtered dataset
  n_filtered <- nrow(filter(dat_filtered, catch_weight_combined>0 | catch_count_combined>0))
  dat_pos2 <- filter(dat_filtered, catch_weight_combined>0 | catch_count_combined>0)
  
  #Make a table
  n_left <- n_total-n_filtered
  n_prop <- n_filtered/n_total
  
  if(i==1){
    dat_lats <- data.frame(common_name=this_species,
                           n_total=n_total,
                           n_filtered=n_filtered,
                           n_left=n_left,
                           n_prop=n_prop)
  } else {
    dat_lats2 <- data.frame(common_name=this_species,
                            n_total=n_total,
                            n_filtered=n_filtered,
                            n_left=n_left,
                            n_prop=n_prop)
    dat_lats <- bind_rows(dat_lats, dat_lats2)
  }
  
  #Disjoin
  temp <- anti_join(dat_pos, dat_pos2)
  
  #Plot
  ggplot(us_coast_proj) + geom_sf() +
  geom_point(temp, mapping=aes(x=X*1000, y=Y*1000), colour="#0D0887FF", size=1)+
  xlim(min(dat2$X)*1000, max(dat2$X)*1000)+
  ylim(min(dat2$Y)*1000, max(dat2$Y)*1000)
  
    
  ggsave(
    paste0("output/", output_folder, "/plots/data_maps/filtered_", this_species, ".png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height = 6,
    units = c("in"),
    dpi = 600,
    limitsize = TRUE, bg="white"
  )
}

#Manually remove limit where do not have clear southern limit
no_s_limits <- c('longspine thornyhead', 'shortspine thornyhead','slender sole', 'spotted ratfish', 'southern rock sole', 'blackbelly eelpout','dover sole', 'english sole', 'lingcod', 'longnose skate', 'pacific hake', 'pacific sanddab', 'petrale sole', 'sablefish')
species_table$southern_limit <- ifelse(species_table$common_name %in% no_s_limits, NA, species_table$southern_limit)

#Manually remove limit where do not have clear Bering Sea limit
no_n_limits <- c("longspine thornyhead", "walleye pollock", "arrowtooth flounder", "flathead sole", "pacific cod", "pacific halibut", "bering skate")
species_table$ebs_limit <- ifelse(species_table$common_name %in% no_n_limits, NA, species_table$ebs_limit)

#Manually edit limit where do not have clear GOA longitudinal limit
species_table$goa_limit <- ifelse(species_table$common_name=="spiny dogfish", -170, species_table$goa_limit)

#Re-run plots and how any removed to check
#Check constraints--see how many positive catches are removed, plot, make sure reasonable
for(i in 1:nrow(species_table)){
  #Filter data for species
  this_species <- species_table$common_name[i]
  dat2 <- filter(dat, common_name==this_species)
  dat_pos <- filter(dat2, catch_weight_combined>0 | catch_count_combined>0)
  #Number of positive catches in complete dataset
  n_total <- nrow(filter(dat2, catch_weight_combined>0 | catch_count_combined>0))
  
  #Filter data to geographic coordinates
  dat_filtered <- dat2
  if(!is.na(species_table$southern_limit[i])){
    dat_filtered <- filter(dat_filtered, latitude >= species_table$southern_limit[i])
  }
  if(!is.na(species_table$ebs_limit[i])){
    #Separate out the Bering Sea only--because Gulf of Alaska messes this up
    dat_filtered_ebs <- filter(dat_filtered, region=="ebs")
    dat_filtered <- filter(dat_filtered, region!="ebs")
    dat_filtered_ebs <- filter(dat_filtered_ebs, latitude <= species_table$ebs_limit[i])
    dat_filtered <- bind_rows(dat_filtered, dat_filtered_ebs)
  }
  if(!is.na(species_table$goa_limit[i])){
    dat_filtered <- filter(dat_filtered, longitude >= species_table$goa_limit[i])
  }
  
  #save data
  saveRDS(dat_filtered, paste0("data/processed_data/fish_filtered/dat_filtered_", this_species, ".rds"))

  #Number of positive catches in filtered dataset
  n_filtered <- nrow(filter(dat_filtered, catch_weight_combined>0 | catch_count_combined>0))
  dat_pos2 <- filter(dat_filtered, catch_weight_combined>0 | catch_count_combined>0)
  
  #Make a table
  n_left <- n_total-n_filtered
  n_prop <- n_filtered/n_total
  
  if(i==1){
    dat_lats <- data.frame(common_name=this_species,
                           n_total=n_total,
                           n_filtered=n_filtered,
                           n_left=n_left,
                           n_prop=n_prop)
  } else {
    dat_lats2 <- data.frame(common_name=this_species,
                            n_total=n_total,
                            n_filtered=n_filtered,
                            n_left=n_left,
                            n_prop=n_prop)
    dat_lats <- bind_rows(dat_lats, dat_lats2)
  }
  
  #Disjoin
  temp <- anti_join(dat_pos, dat_pos2)
  
  #Plot
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(temp, mapping=aes(x=X*1000, y=Y*1000), colour="#0D0887FF", size=1)+
    xlim(min(dat2$X)*1000, max(dat2$X)*1000)+
    ylim(min(dat2$Y)*1000, max(dat2$Y)*1000)
  
  ggsave(
    paste0("output/", output_folder, "/plots/data_maps/filtered_", this_species, ".png"),
    plot = last_plot(),
    device = NULL,
    path = NULL,
    scale = 1,
    width = 8.5,
    height = 6,
    units = c("in"),
    dpi = 600,
    limitsize = TRUE, bg="white"
  )
}

#Save species table with latitudinal limits
write_xlsx(species_table, "data/species_table.xlsx")

#Species with EBS limits
species_ebs <- filter(species_table, !is.na(ebs_limit))$common_name
species_goa <- filter(species_table, !is.na(goa_limit))$common_name
species_cc <- filter(species_table, !is.na(southern_limit))$common_name

#Plots
ggplot(filter(dat_limits, common_name %in% species_ebs & type=="ebs_limit"), aes(x=latitude, y=prop_cumsum_catch))+
  ggtitle("A  Northern Latitudinal Limits")+
  geom_line()+
  geom_vline(filter(species_table, common_name %in% species_ebs), mapping=aes(xintercept=ebs_limit), linetype="dashed")+
  facet_wrap("common_name", labeller=labeller(common_name=label_wrap_gen(15)))+
  ylab("Cumulative Sum of Catch")+
  xlab("Latitude")

ggsave(
  paste0("output/", output_folder, "/plots/cumulative_latitude_combined_N.png"),
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

ggplot(filter(dat_limits, common_name %in% species_goa & type=="goa_limit"), aes(x=longitude, y=prop_cumsum_catch))+
  ggtitle("B  Gulf of Alaska Longitudinal Limits")+
  geom_line()+
  geom_vline(filter(species_table, common_name %in% species_goa), mapping=aes(xintercept=goa_limit), linetype="dashed")+
  facet_wrap("common_name", labeller=labeller(common_name=label_wrap_gen(15)))+
  scale_x_reverse()+
  ylab("Cumulative Sum of Catch")+
  xlab("Latitude")

ggsave(
  paste0("output/", output_folder, "/plots/cumulative_latitude_goa.png"),
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

ggplot(filter(dat_limits, common_name %in% species_cc & type=="southern_limit"), aes(x=latitude, y=prop_cumsum_catch))+
  ggtitle("C Southern Range Limits")+
  geom_line()+
  geom_vline(filter(species_table, common_name %in% species_cc), mapping=aes(xintercept=southern_limit), linetype="dashed")+
  facet_wrap("common_name", labeller=labeller(common_name=label_wrap_gen(15)))+
  scale_x_reverse()+
  ylab("Cumulative Sum of Catch")+
  xlab("Latitude")

ggsave(
  paste0("output/", output_folder, "/plots/cumulative_latitude_cc.png"),
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

#Map data
#Read in all data
files <- list.files(path = "data/processed_data/fish_filtered", pattern = ".rds", full.names=T)
dat <- map(files,readRDS)
dat <- bind_rows(dat)

ggplot(us_coast_proj) + geom_sf() +
  geom_point(dat, mapping=aes(x=X*1000, y=Y*1000), colour="#0D0887FF", size=0.1)+
  ylim(min(dat$Y)*1000, max(dat$Y)*1000)+
  xlab("Longitude")+
  ylab("Latitude")+
  scale_x_continuous(breaks=c(-150,-120), limits=c(min(dat$X)*1000, max(dat$X)*1000))+
  facet_wrap("common_name", ncol=8, labeller=labeller(common_name=label_wrap_gen(8)))+
  theme(panel.spacing.x = unit(1.5, "lines"))

ggsave(
  paste0("output/", output_folder, "/plots/data_by_species.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 12,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)
