remove.packages("sdmTMB")
remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE,  ref="newbreakpt")
library(sdmTMB)
library(sp)
library(dplyr)
library(tidyr)
library(pals)
library(purrr)
library(ggplot2)
library(ggstance)
library(readxl)
library(viridis)
library(ggh4x)
library(ggnewscale)
library(sf)
library(mapview)
library(openxlsx)
library(parallel)
library(ggpubr)
library(stats)
library(patchwork)

set.seed(9881)

setwd("~/Dropbox/GitHub/wsg-choke-species")

#ggplot themes
theme_set(theme_bw(base_size = 16))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

#Load functions
source("code/helper_funs.R")

#List of regions, species
dat_names <- c("cc", "bc", "goa", "ebs", "coastwide", "cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc")
species <- read_excel("data/species_table.xlsx")
species$common_name <- tolower(species$common_name)
species$scientific_name <- tolower(species$scientific_name)
species_table <- species
species <- unique(species$common_name)
species <- unique(species_table$common_name)
species_iphc <- c("sablefish", "pacific cod", "pacific halibut", "yelloweye rockfish", "longnose skate", "big skate", "spiny dogfish", "rougheye rockfish")

#Output folder
output_folder <- "region_comp4"

#Metabolic index models to use
mi_models2use <- c("model13", "model14", "model15")

#All model names
model_names <- c( "model7", "model8", "model13", "model14",  "model15")

#MI Pars
mi_pars <- read.csv("data/taxa_table.csv")
mi_pars$Group <- tolower(mi_pars$Group)

#Taxa lookup
taxa <- read_excel("data/species_table.xlsx")
taxa$MI_Taxa <- tolower(taxa$MI_Taxa)
taxa$common_name <- tolower(taxa$common_name)

#Load aic table for model output
#aic <- as.data.frame(read_excel(paste0("output/",output_folder, "/aic_table.xlsx")))
aic <- readRDS(file = paste0("output/", output_folder, "/aic_table.rds"))
#Make first 6 columns in aic numeric
aic[,1:(ncol(aic)-5)] <- sapply(aic[,1:(ncol(aic)-5)], as.numeric)

##Set up mapping
map_data <- rnaturalearth::ne_countries(scale = "large",
                                        returnclass = "sf",
                                        continent = "North America")

us_coast_proj <- sf::st_transform(map_data, crs = 32610)

###Conditional effect plot
#Manually calculating
##Create dataset used by all
#Load all fish data
files <- list.files(path = "data/processed_data/fish_filtered", pattern = ".rds", full.names=T)
dat <- map(files,readRDS)
dat <- bind_rows(dat)

#Remove any rows with necessary data missing
dat <- dat %>%
  drop_na(depth, mi1, temperature_C, salinity_psu, X, Y, year)

#Remove weird depths
dat <- filter(dat, depth>0)

#Remove oxygen outliers
dat <- filter(dat, O2_umolkg<1500)

##Make prediction grid for conditional effects
#Sequence of oxygen values
dat_pred <- as.data.frame(seq(0,30, length.out=100))
colnames(dat_pred) <- "po2"

#Add other columns for conditional effect fitting
dat_pred$temperature_C <- 12
#Calc invtemp
kelvin = 273.15
boltz = 0.000086173324
tref <- 12
dat_pred$invtemp <- (1 / boltz)  * ( 1 / (tref + 273.15) - 1 / (tref + 273.15))
dat_pred$pred_id <- 1:100

#----Conditional Effect: Manually w/ Monte Carlo parameter draws
for(i in 1:length(species)){
  this_species <- species[i]
  print(this_species)
  preds2 <- dat_pred
  #Calculate metabolic index from correct taxa parameters
  #Pull correct 
  mi_pars <- read.csv("data/taxa_table.csv")
  mi_pars$Group <- tolower(mi_pars$Group)
  this_taxa <- taxa$MI_Taxa[taxa$common_name==this_species]
  mi_pars <- filter(mi_pars,Group==this_taxa)
  
  Eo <- c(mi_pars$Eolow, mi_pars$Eo, mi_pars$Eohigh)
  
  #Calculate metabolic index for species
  preds2$mi1 <- calc_mi(po2=preds2$po2, inv.temp=preds2$invtemp, Eo=Eo[1],fancy=F)
  preds2$mi2 <- calc_mi(po2=preds2$po2, inv.temp=preds2$invtemp, Eo=Eo[2],fancy=F)
  preds2$mi3 <- calc_mi(po2=preds2$po2, inv.temp=preds2$invtemp, Eo=Eo[3],fancy=F)
  
  #List of data types to pull
  this_aic <- filter(aic, species==this_species)
  datnames <- unique(this_aic$`data type`)
  #For each data type, pull all of the models that did fit
  for (j in 1:length(datnames)){
    preds <- preds2
    this_dat <- datnames[j]
    print(this_dat)
    #Pull the data file for scaling
    this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
    #Mean and sd for scaling temp
    mean_temp <- mean(this_datframe$temperature_C, na.rm=T)
    sd_temp <- sd(this_datframe$temperature_C, na.rm=T)
    #Mean and sd for scaling MIs
    mean_mi1 <- mean(this_datframe$mi1, na.rm=T)
    sd_mi1 <- sd(this_datframe$mi1, na.rm=T)
    mean_mi2 <- mean(this_datframe$mi2, na.rm=T)
    sd_mi2 <- sd(this_datframe$mi2, na.rm=T)
    mean_mi3 <- mean(this_datframe$mi3, na.rm=T)
    sd_mi3 <- sd(this_datframe$mi3, na.rm=T)
    #Scale prediction data
    preds$temp_scaled <- (preds$temperature_C-mean_temp)/sd_temp
    preds$temp_scaled2 <- (preds$temp_scaled)^2
    preds$mi1_s <- (preds$mi1-mean_mi1)/sd_mi1
    preds$mi2_s <- (preds$mi2-mean_mi2)/sd_mi2
    preds$mi3_s <- (preds$mi3-mean_mi3)/sd_mi3
    
    #Filter aic table to the datatype to figure out which models to pull
    this_aic <- filter(aic, species==this_species, `data type`==this_dat)
    
    #Just the first 5 columns
    this_aic <- select(this_aic, model_names)
    #Only the columns that are not NAs (which means they didn't pass sanity checks)
    #Remove columns that are NAs
    this_aic <- this_aic %>%
      select(where(~ !any(is.na(.))))
    
    #Get list of these columns (these are the model fits to pull for model averaging)
    models <- colnames(this_aic)
    if(length(models)>0){
      #Pull the model output files
      model_fits <- list()
      for(h in 1:length(models)){
        fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", models[h], ".rds")))
        model_fits[[h]] <- fit
      }
        #Get model weights
        aics <- list()
        for (k in 1:length(models)){
          aic_models <- stats::AIC(model_fits[[k]])
          aics[[k]] <- aic_models
        }
        aics <- unlist(aics)
        delta_aic <- aics - min(aics)
        weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
        set.seed(459384) # for reproducibility and consistency
        for(g in 1:length(model_fits)){
          if(models[g] %in% mi_models2use){
            fit <- model_fits[[g]]
            # Get variance-covariance matrix
            V <- fit[["sd_report"]][["cov.fixed"]]
            #Pull out just threshold
            vcov <- V[grep("threshold", rownames(V)), grep("threshold", colnames(V))]
            slope <- fit[["sd_report"]][["value"]][["s_slope"]]
            thresh <- fit[["sd_report"]][["value"]][["s_cut"]]
          
            #Pull 100 par combinations from variance-covariance matrix
            bp_pars <- MASS::mvrnorm(n=100, mu=c(log(slope),thresh), Sigma=vcov)
            slope_pars <- exp(bp_pars[,1])
            thresh_pars <- bp_pars[,2]
            
            #Calculate the breakpoint effect
            if(models[g]=="model13"){
              mi_s <- preds$mi1_s
            }
            if(models[g]=="model14"){
              mi_s <- preds$mi2_s
            }
            if(models[g]=="model15"){
              mi_s <- preds$mi3_s
            }
            for(l in 1:nrow(bp_pars)){
            pred <- as.data.frame(sapply(mi_s,breakpoint_calc, slope_pars[l],thresh_pars[l]))
            colnames(pred) <- "pred"
            pred$weight <- weights[g]
            pred$sim <- l
            pred$pred_id <-1:100
            #Combine together
            if(l==1){
              pred_model <- pred
            }
            if(l>1){
              pred_model <- bind_rows(pred_model, pred)
            }
            }
      
          } else {
            #For models without pO2 term, conditional effect will be 0--add 100 iterations of this for each prediction point
            pred_model <- as.data.frame(matrix(-1000000, nrow=length(preds$po2)*100))
            colnames(pred_model) <- "pred"
            pred_model$weight <- weights[g]
            pred_model$pred_id <- rep(1:100, 100)
            pred_model$sim <- rep(1:100, each=100)
          }
          if(g==1){
            pred_all <- pred_model
          }
          if(g>1){
            pred_all <- bind_rows(pred_all, pred_model)
          }
        }     
      #Calculate weighted average
      ens_preds <- pred_all %>% 
        group_by(pred_id, sim) %>% 
        summarise(weighted_mean=weighted.mean(exp(pred),weight))  %>% 
        ungroup()
      
      #Ensemble w/ SD
      ens_preds2 <- ens_preds %>% 
        group_by(pred_id) %>% 
        summarise(ensemble_mean=mean(weighted_mean, na.rm=T), ensemble_sd=sd(weighted_mean, na.rm=T)) %>%
        ungroup()
      ens_preds2 <- ens_preds2 %>% 
        mutate(ensemble_lower=(ensemble_mean-ensemble_sd), ensemble_upper=(ensemble_mean+ensemble_sd))%>% 
        ungroup()
      
      #Scaled
      ens_preds3 <- ens_preds2 %>%
        mutate(ensemble_mean_sc= (ensemble_mean/max(ensemble_upper, na.rm=T)),
               ensemble_mean_lower_sc = (ensemble_lower)/max(ensemble_upper, na.rm=T),
               ensemble_mean_upper_sc = (ensemble_upper)/max(ensemble_upper, na.rm=T))
      #Add prediction dataframe back
      preds5 <- left_join(preds, ens_preds3, by="pred_id")
      preds5$species <- this_species
      preds5$data <- this_dat
      if(!grepl("iphc", this_dat)){
      preds5$region <- this_dat
      preds5$data_type <- "bottom trawl only"
      }
      if(grepl("iphc", this_dat)){
      preds5$region <- gsub(" _iphc", "", this_dat)
      preds5$data_type <- "bottom trawl & IPHC"
      }
      if(j==1){
        all_preds <- preds5}
      if(j>1){
        all_preds <- bind_rows(all_preds, preds5)
      }
    }
  }
      
  if(i==1){
    all_preds2 <- all_preds}
  if(i>1){
    all_preds2 <- bind_rows(all_preds2, all_preds)
  }
  }
  
#Save all_preds2
saveRDS(all_preds2, file = paste0("output/", output_folder, "/conditional_effects_data_all.rds"))

##Plot all with SE
#Set region
all_preds2$region <- factor(all_preds2$region, levels=c("ebs", "goa", "bc", "cc", "coastwide"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current", "Coastwide")
names(labs) <- c("ebs", "goa", "bc", "cc", "coastwide")

ggplot(all_preds2, aes(x=po2, y=ensemble_mean_sc))+
  geom_line(aes(colour=region, linetype=data_type))+
 geom_ribbon(aes(ymin=ensemble_mean_lower_sc, ymax=ensemble_mean_upper_sc, fill=region), alpha=0.2)+
  facet_wrap("species", ncol=4, scales="free_y", labeller=labeller(species=label_wrap_gen(15)))+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE, guide="none")+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE, labels=labs )+
  scale_linetype_manual(values=c("dashed", "solid"))+
  guides(fill = guide_legend(nrow = 3, override.aes = list(alpha = 1)), linetype=guide_legend(nrow=2))+
  xlab("Dissolved Oxygen (kPa)")+
  ylab("Effect on Fish Density")

ggsave(
  paste0("output/", output_folder, "/plots/cond_effect_ensemble_all.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height =11,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#Filtering protocol step 1:  aic to only if a 0 in any of columns named in mi_models2use
aic_full <- aic
aic <- filter(aic, aic[,mi_models2use[1]]==0 | aic[,mi_models2use[2]]==0 | aic[,mi_models2use[3]]==0)
species_o2 <- unique(aic$species)

##Filtering protocol step 2: check reasonableness of breakpoint estimates
#Pull 100 draws of parameter estimates (slope and breakpoint) and calculate ensemble
for(i in 1:length(species_o2)){
  this_species <- species_o2[i]
  print(this_species)
  #Calculate metabolic index from correct taxa parameters
  #Pull correct MI values
  mi_pars <- read.csv("data/taxa_table.csv")
  mi_pars$Group <- tolower(mi_pars$Group)
  this_taxa <- taxa$MI_Taxa[taxa$common_name==this_species]
  mi_pars <- filter(mi_pars,Group==this_taxa)
  
  Eo <- c(mi_pars$Eolow, mi_pars$Eo, mi_pars$Eohigh)
  
  #List of data types to pull
  this_aic <- filter(aic, species==this_species)
  datnames <- unique(this_aic$`data type`)
  #For each data type, pull all of the models that did fit
  for (j in 1:length(datnames)){
    this_dat <- datnames[j]
    print(this_dat)
    #Pull the data file for scaling
    this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
    #Mean and sd for unscaling MIs
    mean_mi1 <- mean(this_datframe$mi1, na.rm=T)
    sd_mi1 <- sd(this_datframe$mi1, na.rm=T)
    mean_mi2 <- mean(this_datframe$mi2, na.rm=T)
    sd_mi2 <- sd(this_datframe$mi2, na.rm=T)
    mean_mi3 <- mean(this_datframe$mi3, na.rm=T)
    sd_mi3 <- sd(this_datframe$mi3, na.rm=T)
    
    #Filter aic table to the datatype to figure out which models to pull
    this_aic <- filter(aic, species==this_species, `data type`==this_dat)
    #Just the first 5 columns
    this_aic <- select(this_aic, model_names)
    #Only the columns that are not NAs (which means they didn't pass sanity checks)
    #Remove columns that are NAs
    this_aic <- this_aic %>%
      select(where(~ !any(is.na(.))))
    
    #Get list of these columns (these are the model fits to pull for model averaging)
    models <- colnames(this_aic)
    #Filter to just the MI models
    models <- models[models %in% mi_models2use]
    if(length(models)>0){
      #Pull the model output files
      model_fits <- list()
      for(h in 1:length(models)){
        fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", models[h], ".rds")))
        model_fits[[h]] <- fit
      }
      #Get model weights
      aics <- list()
      for (k in 1:length(models)){
        aic_models <- AIC(model_fits[[k]])
        aics[[k]] <- aic_models
      }
      aics <- unlist(aics)
      delta_aic <- aics - min(aics)
      weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
      set.seed(459384) # for reproducibility and consistency
      for(g in 1:length(model_fits)){
       # if(models[g] %in% mi_models2use){
          fit <- model_fits[[g]]
          pars <- as.data.frame(tidy(fit, effects="fixed", conf.int=T))
          # Get variance-covariance matrix
          V <- fit[["sd_report"]][["cov.fixed"]]
          #Pull out just threshold
          vcov <- V[grep("threshold", rownames(V)), grep("threshold", colnames(V))]
          slope <- fit[["sd_report"]][["value"]][["s_slope"]]
          thresh <- fit[["sd_report"]][["value"]][["s_cut"]]
          #Pull 100 par combinations from variance-covariance matrix
          bp_pars <- MASS::mvrnorm(n=100, mu=c(log(slope),thresh), Sigma=vcov)
          slope_est <- exp(bp_pars[,1])
          bp_par <- bp_pars[,2]
          pars <- data.frame(bp_par=bp_par, slope_est=slope_est)

          if(models[g]=="model13"){
           mean_mi <- mean_mi1
           sd_mi <- sd_mi1
          }
          if(models[g]=="model14"){
            mean_mi <- mean_mi2
            sd_mi <- sd_mi2
          }
          if(models[g]=="model15"){
            mean_mi <- mean_mi3
            sd_mi <- sd_mi3
          }
          pars$bp_pars <- (pars$bp_par*sd_mi)+mean_mi
          pars$sim <- 1:100
          pars$weight <- weights[g]
        # } else {
        #   #For models without pO2 term, conditional effect will be 0--add 100 iterations of this for each prediction point
        #   pars <- as.data.frame(matrix(0, nrow=100))
        #   colnames(pars) <- "bp_pars"
        #   pars$slope_pars <- 0
        #   pars$weight <- weights[g]
        #   pars$sim <- 1:100
       # }
        if(g==1){
          pars_all <- pars
        }
        if(g>1){
          pars_all <- bind_rows(pars_all, pars)
        }
      }
     # }     
      #Calculate weighted average
      pars_weighted <- pars_all %>% 
        group_by(sim) %>% 
        summarise(weighted_mean_bp=weighted.mean(bp_pars,weight, na.rm=T), weighted_mean_slope=weighted.mean(slope_est, weight, na.rm=T))  %>% 
        ungroup()
      
      #Ensemble w/ SD
      pars_ensemble <- pars_weighted %>% 
        summarise(bp_ensemble_mean=mean(weighted_mean_bp, na.rm=T), bp_ensemble_sd=sd(weighted_mean_bp, na.rm=T), slope_ensemble_mean=mean(weighted_mean_slope, na.rm=T), slope_ensemble_sd=sd(weighted_mean_slope, na.rm=T) ) %>%
        ungroup()
      pars_ensemble <- pars_ensemble %>% 
        mutate(ensemble_lower=(bp_ensemble_mean-bp_ensemble_sd), ensemble_upper=(bp_ensemble_mean+bp_ensemble_sd), slope_lower=(slope_ensemble_mean-slope_ensemble_sd), slope_upper=(slope_ensemble_mean+slope_ensemble_sd))%>% 
        ungroup()
      
      #Add other info
      pars_ensemble$species <- this_species
      pars_ensemble$data <- this_dat
      pars_ensemble$max_mi3 <- max(this_datframe$mi3, na.rm=T)
      pars_ensemble$max_mi2 <- max(this_datframe$mi2, na.rm=T)
      pars_ensemble$max_mi1 <- max(this_datframe$mi1, na.rm=T)
      
      if(!grepl("iphc", this_dat)){
        pars_ensemble$region <- this_dat
        pars_ensemble$data_type <- "bottom trawl only"
      }
      if(grepl("iphc", this_dat)){
        pars_ensemble$region <- gsub(" _iphc", "", this_dat)
        pars_ensemble$data_type <- "bottom trawl & IPHC"
      }
      if(j==1){
        all_preds <- pars_ensemble}
      if(j>1){
        all_preds <- bind_rows(all_preds, pars_ensemble)
      }
    }
  }

  if(i==1){
    bp_est <- all_preds}
  if(i>1){
    bp_est<- bind_rows(bp_est, all_preds)
  }
}

##Filter out species that don't have reasonable estimates
bp_est <- unique(bp_est)
bp_est <- filter(bp_est, (bp_ensemble_mean<max_mi1 & ensemble_lower>0 & bp_ensemble_sd<10 & slope_lower>0))

#Creat ID
bp_est$id <- paste(bp_est$species, bp_est$data, sep="_")
all_preds2$id <- paste(all_preds2$species, all_preds2$data, sep="_")

#conditional effects not negative or NA standard errors
all_preds3 <- all_preds2 %>% filter(id %in% bp_est$id)
ids2use <- all_preds3 %>% 
  group_by(id) %>%
  summarize(min=min(ensemble_lower))  %>%
  filter(!is.na(min)&min>0)
bp_est <- filter(bp_est, id %in% ids2use$id)

#Save
saveRDS(bp_est, file = paste0("output/", output_folder, "/breakpoint_estimates.rds"))

#Save a streamlined csv for a table
bp_save <- select(bp_est, species, region, data_type, bp_ensemble_mean, bp_ensemble_sd, slope_ensemble_mean,slope_ensemble_sd)
write.csv(bp_save, file = paste0("output/", output_folder, "/breakpoint_estimates.csv"), row.names=F)

#Combine with aic and save
aic_full$data_type <- ifelse(grepl("iphc", aic_full$`data type`), "bottom trawl & IPHC", "bottom trawl only")
aic_full$region <- case_when(
  grepl("ebs", aic_full$`data type`) ~ "ebs",
  grepl("goa", aic_full$`data type`) ~ "goa",
  grepl("bc", aic_full$`data type`) ~ "bc",
  grepl("cc", aic_full$`data type`) ~ "cc",
  grepl("coastwide", aic_full$`data type`) ~ "coastwide"
)
bp_aic <- left_join(aic_full, bp_save, by=c("species", "data_type", "region"))
bp_aic$"data type" <- NULL
#Round model7-model15, bp_ensemble_mean through slope_ensemble_sd
bp_aic[,c("model7", "model8", "model13", "model14", "model15", "bp_ensemble_mean", "bp_ensemble_sd", "slope_ensemble_mean","slope_ensemble_sd")] <- round(bp_aic[,c("model7", "model8", "model13", "model14", "model15", "bp_ensemble_mean", "bp_ensemble_sd", "slope_ensemble_mean","slope_ensemble_sd")], 3)
#Replace NA with --
bp_aic[is.na(bp_aic)] <- "--"
bp_aic$N_obs <- bp_aic$'N obs'
bp_aic$N_region <- bp_aic$'N region'
bp_aic$N_years <- bp_aic$'N years'
#Reorder
bp_aic <- bp_aic %>%
  select(species, region, data_type, model7, model8, model13, model14, model15, bp_ensemble_mean, bp_ensemble_sd, slope_ensemble_mean, slope_ensemble_sd, N_obs, N_region, N_years)

#save as excel file
write.xlsx(bp_aic, file = paste0("output/", output_folder, "/breakpoint_estimates_aic.xlsx"), rowNames=F)

#Plot, filtered out
ggplot(filter(all_preds2, id %in% bp_est$id), aes(x=po2, y=ensemble_mean_sc))+
  geom_line(aes(colour=region, linetype=data_type))+
  geom_ribbon(aes(ymin=ensemble_mean_lower_sc, ymax=ensemble_mean_upper_sc, fill=region), alpha=0.2)+
  facet_wrap("species", ncol=4, scales="free_y", labeller=labeller(species=label_wrap_gen(15)))+
  theme(legend.title=element_blank(), 
        legend.box.spacing = unit(0, "pt"),
        legend.position="top")+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE, labels=labs)+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE, labels=labs)+
  scale_linetype_manual(values=c("dashed", "solid"))+
  guides(fill = guide_legend(nrow = 2, labels=labs), color=guide_legend(nrow=2, labels=labs, override.aes=list(size=4)), linetype=guide_legend(nrow=2))+
  xlab("Oxygen (kPa) at 12 C")+
  ylab("Effect on Fish Density")

ggsave(
  paste0("output/", output_folder, "/plots/species_filtered_alldat.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 9,
  height= 10,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#Plot line plot
#Set region
bp_est$region <- factor(bp_est$region, levels=c("ebs", "goa", "bc", "cc", "coastwide"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current", "Coastwide")
names(labs) <- c("ebs", "goa", "bc", "cc", "coastwide")

ggplot(bp_est, aes(y=reorder(species, bp_ensemble_mean, decreasing=T), x=bp_ensemble_mean))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=region, shape=data_type), size=3, position=ggstance::position_dodgev(height=0.4))+
  geom_linerange(aes(xmin = ensemble_lower, xmax = ensemble_upper, colour=region),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12))+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  #xlim(-1,10)+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE)+
  guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
  theme(legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0.25, "lines"), #Minimize legend space
        panel.spacing = unit(5, "lines"))+ #Make more space between species
  #scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab("Temperature-Corrected Oxygen Breakpoint Estimate (kPa)")+
  ylab("")
#geom_vline(xintercept=0, linetype="dashed")+
# geom_vline(xintercept=1, linetype="dashed")

ggsave(
  paste0("output/", output_folder, "/plots/breakpoint_estimates_all.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height =8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#All regions, but not all data types
for(i in 1:length(unique(bp_est$species))){
  this_dat <- filter(bp_est, species==unique(bp_est$species)[i])
  this_species <- unique(bp_est$species)[i]
    if(this_species %in% species_iphc){
      this_dat <- filter(this_dat, grepl("iphc", data))
    } 
  if(i==1){
    bp_est_all <- this_dat
  } else {
    bp_est_all <- bind_rows(bp_est_all, this_dat)
  }
}

ggplot(bp_est_all, aes(y=reorder(species, bp_ensemble_mean, decreasing=T), x=bp_ensemble_mean))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=region), size=3, position=ggstance::position_dodgev(height=0.4))+
  geom_linerange(aes(xmin = ensemble_lower, xmax = ensemble_upper, colour=region),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12))+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  #xlim(-1,10)+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE)+
  guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
  theme(legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0.25, "lines"), #Minimize legend space
        panel.spacing = unit(5, "lines"))+ #Make more space between species
  #scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab("Temperature-Corrected Oxygen Breakpoint Estimate (kPa)")+
  ylab("")
#geom_vline(xintercept=0, linetype="dashed")+
# geom_vline(xintercept=1, linetype="dashed")

ggsave(
  paste0("output/", output_folder, "/plots/breakpoint_estimates_regions.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height =8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

##Plot only one data type per species
#IPHC if available, coastwide if available
for(i in 1:length(unique(bp_est$species))){
  this_dat <- filter(bp_est, species==unique(bp_est$species)[i])
  this_species <- unique(bp_est$species)[i]
  if(any(grepl("coastwide", this_dat$data))){
    this_dat <- filter(this_dat, grepl("coastwide", data))
    if(this_species %in% species_iphc){
      this_dat <- filter(this_dat, grepl("iphc", data))
    } 
  }
  if(i==1){
    bp_est2 <- this_dat
  } else {
    bp_est2 <- bind_rows(bp_est2, this_dat)
  }
}

saveRDS(bp_est2, file = paste0("output/", output_folder, "/breakpoint_estimates_filtered.rds"))

#Plot
ggplot(filter(all_preds2, id %in% bp_est2$id), aes(x=po2, y=ensemble_mean_sc))+
  geom_line(aes(colour=region, linetype=data_type))+
  geom_ribbon(aes(ymin=ensemble_mean_lower_sc, ymax=ensemble_mean_upper_sc, fill=region), alpha=0.2)+
  facet_wrap("species", ncol=4, scales="free_y")+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE, labels=labs )+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE,labels=labs )+
  scale_linetype_manual(values=c("dashed", "solid"))+
  guides(fill = guide_legend(nrow = 2, labels=labs), color=guide_legend(nrow=2, labels=labs, override.aes=list(size=4)), linetype=guide_legend(nrow=2))+
  xlab("Oxygen (kPa) at 12 C")+
  ylab("Effect on Fish Density")

ggsave(
  paste0("output/", output_folder, "/plots/species_filtered_onedat.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =9.5,
  height =8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

##Line plot
bp_est$id <- paste(bp_est$species, bp_est$data, sep="_")
ggplot(bp_est2, aes(y=reorder(species, bp_ensemble_mean, decreasing=T), x=bp_ensemble_mean, colour=region))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=region, shape=data_type), size=3, position=ggstance::position_dodgev(height=0.4))+
  #Can add back shape
  geom_linerange(aes(xmin = ensemble_lower, xmax = ensemble_upper, colour=region),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12))+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  #xlim(-1,10)+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE)+
  guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
  theme(legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0.25, "lines"), #Minimize legend space
        panel.spacing = unit(5, "lines"))+ #Make more space between species
  #scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab("Temperature-Corrected Oxygen Threshold (kPa)")+
  ylab("")
#geom_vline(xintercept=0, linetype="dashed")+
# geom_vline(xintercept=1, linetype="dashed")

ggsave(
  paste0("output/", output_folder, "/plots/threshold_est_filtered_onedat.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =9.5,
  height =8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#Just one region per species
for(i in 1:length(unique(bp_est2$species))){
  this_dat <- filter(bp_est2, species==unique(bp_est2$species)[i])
  this_species <- unique(bp_est2$species)[i]
    if(nrow(this_dat)>1){
      this_dat <- filter(this_dat, region=="bc")
    }
  if(i==1){
    bp_est_temp <- this_dat
  } else {
    bp_est_temp <- bind_rows(bp_est_temp, this_dat)
  }
}

ggplot(bp_est_temp, aes(y=reorder(species, bp_ensemble_mean, decreasing=T), x=bp_ensemble_mean, colour=region))+
  #facet_grid(rows="species", scales="free_y", space="free_y", switch="y")+
  geom_point(aes(colour=region), size=3, position=ggstance::position_dodgev(height=0.4))+
  #Can add back shape
  geom_linerange(aes(xmin = ensemble_lower, xmax = ensemble_upper, colour=region),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12))+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=15))+
  #xlim(-1,10)+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677", "#332288"), drop=FALSE)+
  guides( color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))+
  theme(legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0.25, "lines"), #Minimize legend space
        panel.spacing = unit(5, "lines"))+ #Make more space between species
  #scale_shape_discrete(labels=c("low Eo", "median Eo", "high Eo"))+
  xlab("Temperature-Corrected Oxygen Threshold (kPa)")+
  ylab("")
#geom_vline(xintercept=0, linetype="dashed")+
# geom_vline(xintercept=1, linetype="dashed")

ggsave(
  paste0("output/", output_folder, "/plots/threshold_est_filtered_oneregion.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =9.5,
  height =8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

#With only the single selected region/data (in black), and in ascending order
#Keep only california current for regions
ggplot(bp_est_temp, aes(y=reorder(species, bp_ensemble_mean, decreasing=T), x=bp_ensemble_mean))+
  geom_point(size=3, position=ggstance::position_dodgev(height=0.4))+
  geom_linerange(aes(xmin = ensemble_lower, xmax = ensemble_upper),  position=ggstance::position_dodgev(height=0.4), size=1, alpha=0.5)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=12))+
  theme(legend.position="top")+
  theme(legend.title=element_blank())+
  theme(text=element_text(size=20))+
  #xlim(-1,10)+
  theme(legend.box = "vertical",
        legend.spacing.y = unit(0, "pt"),
        legend.key.height = unit(0.25, "lines"), #Minimize legend space
        panel.spacing = unit(5, "lines"))+ #Make more space between species
  xlab("Temperature-Adjusted Oxygen Threshold (kPa)")+
  ylab("")

ggsave(
  paste0("output/", output_folder, "/plots/thresholds_lineplot_for_presentation.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =9.5,
  height =8.5,
  units = c("in"),
  bg="white",
  dpi = 600,
  limitsize = TRUE
)

###Evaluate depth effects
#Depth sequence
for(i in 1:nrow(bp_est)) {
  this_species <- bp_est$species[i]
  this_data <- bp_est$data[i]
  print(this_species)
  print(this_data)
  #Pull model fit
  this_fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_data,"_model15.rds")))
  this_dat <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_data, "_dat.rds")))
  if(class(this_fit)!="try-error"){
    pars <- as.data.frame(tidy(this_fit, effects="fixed", conf.int=T))
    depth1 <- filter(pars,term=="log_depth_scaled")
    depth1 <- depth1$estimate
    depth2 <- filter(pars,term=="log_depth_scaled2")
    depth2 <- depth2$estimate
    depth3 <- filter(pars,term=="log_depth_scaled3")
    depth3 <- depth3$estimate
    
    #Sequence of depths
    depth_seq <- as.data.frame(seq(min(this_dat$log_depth_scaled, na.rm=T), max(this_dat$log_depth_scaled, na.rm=T), length.out=100))
    colnames(depth_seq) <- "depth_seq1"
    depth_seq$depth_seq2 <- depth_seq$depth_seq1^2
    depth_seq$depth_seq3 <- depth_seq$depth_seq1^3
    
    #Calculate conditional effect
    depth_seq$effect <- depth1*depth_seq$depth_seq1 + depth2*depth_seq$depth_seq2 + depth3*depth_seq$depth_seq3
 
    #Plot and save
    ggplot(depth_seq, aes(x=depth_seq1, y=exp(effect)))+
      geom_line()+
      xlab("Log Depth (scaled)")+
      ylab("Conditional Effect on Fish Density")+
      ggtitle(paste(this_species, this_data))+
      theme(text=element_text(size=15))
    
    ggsave(
      paste0("output/", output_folder, "/depth_conditional_effects/", this_species, "_", this_data, "_depth_conditional_effects.png"),
      plot = last_plot(),
      device = NULL,
      path = NULL,
      scale = 1,
      width = 8.5,
      height = 5,
      units = c("in"),
      bg="white",
      dpi = 600,
      limitsize = TRUE
    )
  }
}

##-----Historical observations and conditional effects
bp_est2 <- readRDS(paste0("output/", output_folder, "/breakpoint_estimates_filtered.rds"))

###Loop through for each species for all data (limited to range)
#----Conditional Effect: Manually w/ Monte Carlo parameter draws
for(i in 1:nrow(bp_est2)){
  this_species <- bp_est2$species[i]
  this_data <- bp_est2$data[i]
  print(this_species)
  #Calculate metabolic index from correct taxa parameters
  #Pull correct 
  species_tab <- filter(species_table, common_name==this_species)
  
  #Pull the data file
  this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_data, "_dat.rds")))
  
  #Make prediction dataframe
  preds <- this_datframe
  
  #Add pred_id for re-merging with predictions
  preds$pred_id <- 1:nrow(preds)
  
  #Mean and sd for unscaling MIs
  mean_mi1 <- mean(this_datframe$mi1, na.rm=T)
  sd_mi1 <- sd(this_datframe$mi1, na.rm=T)
  mean_mi2 <- mean(this_datframe$mi2, na.rm=T)
  sd_mi2 <- sd(this_datframe$mi2, na.rm=T)
  mean_mi3 <- mean(this_datframe$mi3, na.rm=T)
  sd_mi3 <- sd(this_datframe$mi3, na.rm=T)
  
  #Calculate
  preds$mi1_s <- (preds$mi1-mean_mi1)/sd_mi1
  preds$mi2_s <- (preds$mi2-mean_mi2)/sd_mi2
  preds$mi3_s <- (preds$mi3-mean_mi3)/sd_mi3

  #List of data types to pull
  this_aic <- filter(aic, species==this_species& `data type`==this_data)
  #For each data type, pull all of the models that did fit
  this_dat <- this_data
  print(this_dat)
  #Just the first 5 columns
  this_aic <- select(this_aic, model_names)
  #Only the columns that are not NAs (which means they didn't pass sanity checks)
  #Remove columns that are NAs
  this_aic <- this_aic %>%
    select(where(~ !any(is.na(.))))
    
  #Get list of these columns (these are the model fits to pull for model averaging)
  models <- colnames(this_aic)
  
  #Pull the model output files
  model_fits <- list()
  for(h in 1:length(models)){
    fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", models[h], ".rds")))
    model_fits[[h]] <- fit
  }
  #Get model weights
  aics <- list()
  for (k in 1:length(models)){
  aic_models <- stats::AIC(model_fits[[k]])
  aics[[k]] <- aic_models
  }
 aics <- unlist(aics)
 delta_aic <- aics - min(aics)
 weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
 set.seed(459384) # for reproducibility and consistency
 for(g in 1:length(model_fits)){
        if(models[g] %in% mi_models2use){
          fit <- model_fits[[g]]
          # Get variance-covariance matrix
          V <- fit[["sd_report"]][["cov.fixed"]]
          #Pull out just threshold
          vcov <- V[grep("threshold", rownames(V)), grep("threshold", colnames(V))]
          slope <- fit[["sd_report"]][["value"]][["s_slope"]]
          thresh <- fit[["sd_report"]][["value"]][["s_cut"]]
          
          #Pull 100 par combinations from variance-covariance matrix
          bp_pars <- MASS::mvrnorm(n=100, mu=c(log(slope),thresh), Sigma=vcov)
          slope_pars <- exp(bp_pars[,1])
          thresh_pars <- bp_pars[,2]
          
          #Calculate the breakpoint effect
          if(models[g]=="model13"){
            mi_s <- preds$mi1_s
          }
          if(models[g]=="model14"){
            mi_s <- preds$mi2_s
          }
          if(models[g]=="model15"){
            mi_s <- preds$mi3_s
          }
          for(l in 1:nrow(bp_pars)){
            pred <- as.data.frame(sapply(mi_s,breakpoint_calc, slope_pars[l],thresh_pars[l]))
            colnames(pred) <- "pred"
            #repeat this value the length of nrow(dat)
            max <- rep((thresh_pars[l]+1), nrow(preds))
            pred$max_effect <- sapply(max,breakpoint_calc, slope_pars[l],thresh_pars[l])
            pred$est_effect_prop <- (exp(pred$max_effect)-exp(pred$pred))/exp(pred$max_effect)
            pred$weight <- weights[g]
            pred$sim <- l
            pred$pred_id <-1:nrow(preds)
            
            #Combine together
            if(l==1){
              pred_model <- pred
            }
            if(l>1){
              pred_model <- bind_rows(pred_model, pred)
            }
          }
          
        } else {
          #For models without pO2' term, conditional effect will be 0--add 100 iterations of this for each prediction point
          pred_model <- as.data.frame(matrix(0, nrow=nrow(preds)*100))
          colnames(pred_model) <- "pred"
          pred_model$est_effect_prop <- 0
          pred_model$weight <- weights[g]
          pred_model$pred_id <- rep(1:nrow(preds), 100)
          pred_model$sim <- rep(1:100, each=nrow(preds))
        }
        if(g==1){
          pred_all <- pred_model
        }
        if(g>1){
          pred_all <- bind_rows(pred_all, pred_model)
        }
 }     
    
  #Calculate weighted average
      ens_preds <- pred_all %>% 
        group_by(pred_id, sim) %>% 
        summarise(weighted_mean=weighted.mean(est_effect_prop,weight, na.rm=T))  %>% 
        ungroup()
      
      #Ensemble w/ SD
      ens_preds2 <- ens_preds %>% 
        group_by(pred_id) %>% 
        summarise(ensemble_mean=mean(weighted_mean, na.rm=T), ensemble_sd=sd(weighted_mean, na.rm=T)) %>%
        ungroup()
      ens_preds2 <- ens_preds2 %>% 
        mutate(ensemble_lower=(ensemble_mean-ensemble_sd), ensemble_upper=(ensemble_mean+ensemble_sd))%>% 
        ungroup()
      
      #Add prediction dataframe back
      preds5 <- left_join(preds, ens_preds2, by="pred_id")
      #Other columns
      preds5$data <- this_data 
      preds5$common_name <- this_species
      preds5$depth_limit <- species_tab$depth
      preds5$min_lat <- species_tab$northern_limit
      preds5$max_lat <- species_tab$southern_limit
      if(!grepl("iphc", this_dat)){
        preds5$data_type <- "bottom trawl only"
      }
      if(grepl("iphc", this_dat)){
        preds5$data_type <- "bottom trawl & IPHC"
      }
      if(i==1){
        all_obs <- preds5}
      if(i>1){
        all_obs <- bind_rows(all_obs, preds5)
      }
}

saveRDS(all_obs, file=paste0("output/", output_folder, "/", "conditional_effects_obs_points_all.rds"))

##Plot just coastwide species, on one plot
dat2plot <- filter(all_obs, grepl("coastwide", data))

ggplot(us_coast_proj) + geom_sf() +
  #geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="#0D0887FF", size=0.5, alpha=0.1)+
  geom_point(filter(dat2plot, ensemble_mean==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.1, alpha=0.1)+
  geom_point(filter(dat2plot,ensemble_mean>0),mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.1, alpha=0.1)+
   #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  facet_wrap("common_name", labeller=labeller(common_name=label_wrap_gen(20)))+
  scale_x_continuous(breaks=c(-150,-135,-120), limits=c(min(dat2plot$X)*1000, max(dat2plot$X)*1000))+
  ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position=c(0.8,0.15), legend.justification="center",
        legend.box.spacing = unit(0, "pt"), panel.spacing = unit(1.2, "lines"))+
  scale_colour_viridis(name="Reduction in \n local density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/map_coastwide_species_combined.png"),
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

###Add full predictions
###Add full density predictions and density predictions with O2 above threshold
coastwide_full <- list()
species_list <- unique(dat2plot$common_name)
#Calculate predictions for each grid cell
for(i in 1:length(species_list)){
  this_species <- species_list[i]
  bp_est3 <- filter(bp_est2, species==this_species)
  this_data <- bp_est3$data
  print(this_species)
  #Calculate metabolic index from correct taxa parameters
  #Pull correct 
  species_tab <- filter(species_table, common_name==this_species)
  #Calculate metabolic index from correct taxa parameters
  this_taxa <- species_tab$MI_Taxa[species_tab$common_name==this_species]
  mi_pars2 <- filter(mi_pars,Group==this_taxa)
  
  Eo <- c(mi_pars2$Eolow, mi_pars2$Eo, mi_pars2$Eohigh)
  
  #New dataframe
  preds <- filter(dat2plot, common_name==this_species)
  
  #Prediction ID (row number)
  preds$pred_id <- 1:nrow(preds)
  
  #Which aic/data to pull
  this_aic <- filter(aic, species==this_species& `data type`==this_data)
  this_dat <- this_data
  print(this_dat)
  #Just the first 5 columns
  this_aic <- select(this_aic, model_names)
  #Only the columns that are not NAs (which means they didn't pass sanity checks)
  #Remove columns that are NAs
  this_aic <- this_aic %>%
    select(where(~ !any(is.na(.))))
  
  #Get list of these columns (these are the model fits to pull for model averaging)
  models <- colnames(this_aic)
  
  #Pull the model output files
  model_fits <- list()
  for(h in 1:length(models)){
    fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", models[h], ".rds")))
    model_fits[[h]] <- fit
  }
  #Get model weights
  aics <- list()
  for (k in 1:length(models)){
    aic_models <- stats::AIC(model_fits[[k]])
    aics[[k]] <- aic_models
  }
  aics <- unlist(aics)
  delta_aic <- aics - min(aics)
  weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
  set.seed(459384) # for reproducibility and consistency
  
  for(g in 1:length(model_fits)){
    #Calculate the breakpoint effect
    #Prediction
    print(models[[g]])
    prediction <- predict(fit, newdata=preds, se.fit=FALSE)
    prediction$weight <- weights[g]
    if(g==1){
      df <- prediction
    }
    if(g>1){
      df <- bind_rows(df, prediction)
    }
  }     
  
  #Calculate weighted average
  ens_preds <- df %>% 
    group_by(pred_id) %>% 
    summarise(weighted_mean=weighted.mean(est,weight, na.rm=T))  %>% 
    ungroup()
  
  #Ensemble w/ SD
  ens_preds2 <- ens_preds %>% 
    group_by(pred_id) %>% 
    summarise(ensemble_mean=mean(weighted_mean, na.rm=T))
  
  #Add prediction dataframe back
  preds5 <- left_join(preds, ens_preds2, by="pred_id")
  
  #And for super high O2
  #Add dummy MI columns for MI above threshold
  #Keep real values
  preds$mi1_s_real <- preds$mi1_s
  preds$mi2_s_real <- preds$mi2_s
  preds$mi3_s_real <- preds$mi3_s
  preds$mi1_s <- 1000
  preds$mi2_s <- 1000
  preds$mi3_s <- 1000
  
  for(g in 1:length(model_fits)){
    #Calculate the breakpoint effect
    #Prediction
    prediction <- predict(fit, newdata=preds, se.fit=FALSE)
    prediction$weight <- weights[g]
    if(g==1){
      df <- prediction
    }
    if(g>1){
      df <- bind_rows(df, prediction)
    }
  }     
  
  #Calculate weighted average
  ens_preds <- df %>% 
    group_by(pred_id) %>% 
    summarise(weighted_mean=weighted.mean(est,weight, na.rm=T))  %>% 
    ungroup()
  
  #Ensemble w/ SD
  ens_preds2 <- ens_preds %>% 
    group_by(pred_id) %>% 
    summarise(ensemble_mean_full=mean(weighted_mean, na.rm=T))
  
  #Add back
  #Add prediction dataframe back
  preds5 <- left_join(preds5, ens_preds2, by="pred_id")
  
  #Scale full density
  preds5$density_scaled <- exp(preds5$ensemble_mean_full)/exp(max(preds5$ensemble_mean_full))
  
  #Other columns
  preds5$data <- this_data 
  preds5$common_name <- this_species
  preds5$depth_limit <- species_tab$depth
  preds5$min_lat <- species_tab$southern_limit
  
  #Make MI back to normal
  preds5$mi1_s <- preds$mi1_s_real
  preds5$mi2_s <- preds$mi2_s_real
  preds5$mi3_s <- preds$mi3_s_real
  #remove fake MI
  preds5$mi1_s_real <- NULL
  preds5$mi2_s_real <- NULL
  preds5$mi3_s_real <- NULL
  
  if(!grepl("iphc", this_dat)){
    preds5$data_type <- "bottom trawl only"
  }
  if(grepl("iphc", this_dat)){
    preds5$data_type <- "bottom trawl & IPHC"
  }
  saveRDS(preds5, file = paste0("output/", output_folder, "/", this_species, "_", "coastwide_preds.rds"))
  coastwide_full[[i]] <- preds5
}

#Plot scatterplot
coastwide_preds <- bind_rows(coastwide_full)
ggplot(coastwide_preds, aes(x=density_scaled, y=ensemble_mean.x))+
  geom_point(alpha=0.1, size=0.01)+
  facet_wrap("common_name", scales="free", labeller=labeller(common_name=label_wrap_gen(20)))+
  xlab("Scaled Predicted Density")+
  ylab("Reduction in Local Density")+
  theme(legend.position="top",
        legend.direction = "horizontal",
        legend.title.align = 0.5,
        legend.title.position = "top",
        legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_x_continuous(labels=c(0,.25,.5,.75,1))

ggsave(
  paste0("output/", output_folder, "/", "plots/coastwide_scatterplot.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 13,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#Colored by depth
ggplot(coastwide_preds, aes(x=density_scaled, y=ensemble_mean.x))+
  geom_point(alpha=0.1, size=0.01, aes(color=depth))+
  facet_wrap("common_name", scales="free", labeller=labeller(common_name=label_wrap_gen(20)))+
  xlab("Scaled Predicted Density")+
  ylab("Reduction in Local Density")+
  theme(legend.position="top",
        legend.direction = "horizontal",
        legend.title.align = 0.5,
        legend.title.position = "top",
        legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_x_continuous(labels=c(0,.25,.5,.75,1))

ggsave(
  paste0("output/", output_folder, "/", "plots/coastwide_scatterplot_depth.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 13,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#Color by latitude
ggplot(coastwide_preds, aes(x=density_scaled, y=ensemble_mean.x))+
  geom_point(alpha=0.1, size=0.01, aes(color=latitude))+
  facet_wrap("common_name", scales="free", labeller=labeller(common_name=label_wrap_gen(20)))+
  xlab("Scaled Predicted Density")+
  ylab("Reduction in Local Density")+
  theme(legend.position="top",
        legend.direction = "horizontal",
        legend.title.align = 0.5,
        legend.title.position = "top",
        legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_x_continuous(labels=c(0,.25,.5,.75,1))

ggsave(
  paste0("output/", output_folder, "/", "plots/coastwide_scatterplot_latitude.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 13,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#Region
ggplot(coastwide_preds, aes(x=density_scaled, y=ensemble_mean.x))+
  geom_point(alpha=0.1, size=0.01, aes(color=region))+
  facet_wrap("common_name", scales="free", labeller=labeller(common_name=label_wrap_gen(20)))+
  xlab("Scaled Predicted Density")+
  ylab("Reduction in Local Density")+
  theme(legend.position="top",
        legend.direction = "horizontal",
        legend.title.align = 0.5,
        legend.title.position = "top",
        legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_x_continuous(labels=c(0,.25,.5,.75,1))

ggsave(
  paste0("output/", output_folder, "/", "plots/coastwide_scatterplot_region.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 13,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Plot scaled biomass
ggplot(us_coast_proj) + geom_sf() +
  #geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="#0D0887FF", size=0.5, alpha=0.1)+
  geom_point(coastwide_preds,mapping=aes(x=X*1000, y=Y*1000, colour=density_scaled), size=0.01, alpha=0.1)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  facet_wrap("common_name", labeller=labeller(common_name=label_wrap_gen(20)))+
  scale_x_continuous(breaks=c(-150,-135,-120), limits=c(min(dat2plot$X)*1000, max(dat2plot$X)*1000))+
  ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position=c(0.8,0.15), legend.justification="center",
        legend.box.spacing = unit(0, "pt"), panel.spacing = unit(1.2, "lines"))+
  scale_colour_viridis(name="Scaled \n local density ", breaks=c(0,0.25, 0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish, option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/map_coastwide_local_density.png"),
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

##Plot coastwide model, but only show points above certain biomass threshold
ggplot(us_coast_proj) + geom_sf() +
  #geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="#0D0887FF", size=0.5, alpha=0.1)+
  geom_point(filter(coastwide_preds, ensemble_mean.x==0&density_scaled<0.1),mapping=aes(x=X*1000, y=Y*1000), colour="grey", alpha=0.1,size=0.01)+
  geom_point(filter(coastwide_preds,ensemble_mean.x>0),mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean.x), size=0.01)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  facet_wrap("common_name", labeller=labeller(common_name=label_wrap_gen(20)))+
  scale_x_continuous(breaks=c(-150,-135,-120), limits=c(min(dat2plot$X)*1000, max(dat2plot$X)*1000))+
  ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position=c(0.8,0.15), legend.justification="center",
        legend.box.spacing = unit(0, "pt"), panel.spacing = unit(1.2, "lines"))+
  scale_colour_viridis(name="Reduction in \n local density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/map_coastwide_combined_alpha_density.png"),
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
##Plot annually
toplot <- unique(dat2plot$common_name)
for(i in 1:length(toplot)) {
  this_species <- toplot[i]
  temp <- filter(dat2plot, common_name==this_species)
  ggplot(us_coast_proj) + geom_sf() +
    #geom_point(filter(dat2plot, est_effect_prop==0),mapping=aes(x=X*1000, y=Y*1000), colour="#0D0887FF", size=0.5, alpha=0.1)+
    geom_point(filter(temp, ensemble_mean==0),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.5, alpha=0.1)+
    geom_point(filter(temp,ensemble_mean>0),mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.9)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    facet_wrap("year", ncol=4)+
    scale_x_continuous(breaks=c(-150,-130), limits=c(min(temp$X)*1000, max(temp$X)*1000))+
    ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
    #facet_wrap("year", ncol=5)+
    theme_minimal(base_size=18)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(legend.position=c(0.92,0.1), panel.spacing = unit(1, "lines"))+
    ggtitle(this_species)+
    scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")
  #  limits=c(0,1), breaks=c(0,0.25, 0.5, 0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish)
  
  ggsave(
    paste0("output/", output_folder, "/", "plots/points_full_annual/map_",this_species,".png"),
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

#Plot observations with example species
dat2use <- filter(all_obs, common_name=="pacific halibut"|common_name=="pacific cod")
dat2use$region <- factor(dat2use$region, levels=c("ebs", "goa", "bc", "cc"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current")
dat2use <- filter(dat2use, year=="2011"|year=="2015"|year=="2023")
names(labs) <- c("ebs", "goa", "bc", "cc")
bp_est2$common_name <- bp_est2$species
ggplot(dat2use, aes(x=po2))+
  geom_density(aes(colour=region))+
  facet_grid(year~common_name)+
  theme(legend.position="top")+
  xlab("Oxygen (kPa)")+
guides(color = guide_legend(nrow = 2))+
  geom_vline(filter(bp_est2, species %in% unique(dat2use$common_name)&(region=="goa"|region=="coastwide")), mapping=aes(xintercept=bp_ensemble_mean), linetype="dashed")+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs, name="")+
  ylab("Density")


ggsave(
  paste0("output/", output_folder, "/", "plots/density_obs_example.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8,
  height = 7,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##proportion of observations in each region with decline
dat2plot <- drop_na(dat2plot, ensemble_mean, ensemble_lower, ensemble_upper)
effects_full_annual_summary <- dat2plot %>%
  group_by(year, region, common_name) %>%
  summarise(prop_below = sum(ensemble_mean>.1)/n(), N=n()) %>%
  ungroup()
effects_full_annual_summary2 <- dat2plot %>%
  group_by(year, region, common_name) %>%
  summarise(prop_below = sum(ensemble_lower>.1)/n()) %>%
  ungroup()
effects_full_annual_summary3 <- dat2plot %>%
  group_by(year, region, common_name) %>%
  summarise(prop_below = sum(ensemble_upper>.1)/n()) %>%
  ungroup()

effects_full_annual_summary$se1 <- effects_full_annual_summary2$prop_below
effects_full_annual_summary$se2 <- effects_full_annual_summary3$prop_below

#Factor region
effects_full_annual_summary$region <- factor(effects_full_annual_summary$region, levels=c("ebs", "goa", "bc", "cc"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current")
names(labs) <- c("ebs", "goa", "bc", "cc")

#At least 50 observations in region
effects_full_annual_summary <- filter(effects_full_annual_summary, N>50)
#plot
ggplot(effects_full_annual_summary, aes(x=year, y=prop_below, group=region))+
  geom_line(aes(colour=region, group=region))+
  #xlim(2010,2024)+
  geom_ribbon(aes(ymin=se2, ymax=se1, fill=region, group=region), alpha=0.2)+
  facet_wrap("common_name", ncol=3)+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=20))+
 scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs )+
 scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  xlab("Year") +
  ylab("Proportion of observations \n w/ >10% local density reduction")

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_below_obs_annual.png"),
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

####Proportion below threshold by depth
dat2plot <- drop_na(dat2plot, ensemble_mean, ensemble_lower, ensemble_upper)
dat2plot$depth_bin <- cut(dat2plot$depth, labels=FALSE, breaks=seq(0, 1500, by=50), include.lowest = TRUE)
dat2plot$depth_bin <- dat2plot$depth_bin*50
effects_full_annual_summary <- dat2plot %>%
  group_by(depth_bin, region, common_name) %>%
  summarise(prop_below = sum(ensemble_mean>0.1)/n(), N=n()) %>%
  ungroup()
effects_full_annual_summary2 <- dat2plot %>%
  group_by(depth_bin, region, common_name) %>%
  summarise(prop_below = sum(ensemble_lower>0.1)/n()) %>%
  ungroup()
effects_full_annual_summary3 <- dat2plot %>%
  group_by(depth_bin, region, common_name) %>%
  summarise(prop_below = sum(ensemble_upper>0.1)/n()) %>%
  ungroup()
effects_full_annual_summary$se1 <- effects_full_annual_summary2$prop_below
effects_full_annual_summary$se2 <- effects_full_annual_summary3$prop_below


#Factor region
effects_full_annual_summary$region <- factor(effects_full_annual_summary$region, levels=c("ebs", "goa", "bc", "cc"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current")
names(labs) <- c("ebs", "goa", "bc", "cc")
effects_full_annual_summary <- filter(effects_full_annual_summary, N>50)
ggplot(effects_full_annual_summary, aes(x=depth_bin, y=prop_below))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", scales="free_x")+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=20))+
  scale_colour_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  xlab("Depth (m)") +
  ylab("Proportion of observations w/ >10% local density reduction")+
  geom_vline(dat2plot, mapping=aes(xintercept=depth_limit), linetype="dashed", colour="black")

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_below_by_depth.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

###Predict to grid for California Current and British Columbia only
bc_grid <- readRDS("data/processed_data/o2/bc_predictions_grid.rds")
bc_grid$latitude <- bc_grid$latitute
cc_grid <- readRDS("data/processed_data/o2/cc_predictions_grid.rds")
cc_grid$latitude <- cc_grid$latitute

#Calculate inverse temp
#Calc invtemp
kelvin = 273.15
boltz = 0.000086173324
tref <- 12
cc_grid$invtemp <- (1 / boltz)  * ( 1 / (cc_grid$temperature_C + 273.15) - 1 / (tref + 273.15))
bc_grid$invtemp <- (1 / boltz)  * ( 1 / (bc_grid$temperature_C + 273.15) - 1 / (tref + 273.15))
cc_grid <- filter(cc_grid, survey=="NWFSC.Combo")
#combine grid
grid <- bind_rows(bc_grid, cc_grid)

#Set mapping extent
ylims <- c(min(grid$Y)*1000, max(grid$Y)*1000)
xlims <- c(min(grid$X)*1000, max(grid$X)*1000)

##Plot oxygen data
#Just WA
grid_wa <- filter(grid, latitude>46&latitude<48.5)
ylims2 <- c(min(grid_wa$Y)*1000, max(grid_wa$Y)*1000)
xlims2 <- c(min(grid_wa$X)*1000, max(grid_wa$X)*1000)
ggplot(us_coast_proj) + geom_sf() +
geom_point(grid_wa,mapping=aes(x=X*1000, y=Y*1000, colour=po2), size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-125), limits=c(xlims2))+
  ylim(ylims2)+
  facet_wrap("year")+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position=c(0.7,0.1), panel.spacing = unit(1.5, "lines"), strip.text=element_text(size=13))+
 scale_colour_viridis(option="mako",limits=c(0,3),oob=scales::squish, name="Oxygen (kPa)", labels=c(0,1,2,">3"))

ggsave(
  paste0("output/", output_folder, "/", "oxygen_grid_wa_annual.png"),
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

options(dplyr.summarise.inform = FALSE)

#Remove goa from spiny dogfish
bp_est2 <- filter(bp_est2, !(species=="spiny dogfish"&data=="goa _iphc"))

#Function to calculate for each grid cell
calc_grid <- function(mi1_s, mi2_s, mi3_s, pred_id, pars_list, models, model_fits, mi_models2use, weights,n) {
  # Define locally within the worker
  breakpoint_calc <- function(x, b_slope, b_thresh){
    if (x < b_thresh) {
      pred = x  * b_slope
      #pred = (x-b_thresh) * b_slope
    } else {
      pred=b_thresh * b_slope
      #pred = 0
    }
    return(pred)
  }
  
  for(g in 1:length(model_fits)){
    if(models[g] %in% mi_models2use){
      #Calculate the breakpoint effect
      if(models[g]=="model13"){
        mi_s <- mi1_s
      }
      if(models[g]=="model14"){
        mi_s <- mi2_s
      }
      if(models[g]=="model15"){
        mi_s <- mi3_s
      }
      #Pull bp_pars and slope_pars from list
      thresh_pars <- pars_list[[g]]$thresh_pars
      slope_pars <- pars_list[[g]]$slope_pars
      
      for(l in 1:length(thresh_pars)){
        pred <- sapply(mi_s,breakpoint_calc, slope_pars[l],thresh_pars[l])
        #repeat this value the length of nrow(dat)
        max <- thresh_pars[l]+1
        max_effect <- sapply(max,breakpoint_calc, slope_pars[l],thresh_pars[l])
        est_effect_prop <- (exp(max_effect)-exp(pred))/exp(max_effect)
        weight <- weights[g]
        sim <- l
        pred_id <- pred_id
        
        #Combine together
        if(l==1){
          est_effect_prop2 <- est_effect_prop
          weight2 <- weight
          sim2 <- sim
          pred_id2 <- pred_id
        }
        
        if(l>1){
          est_effect_prop2 <- c(est_effect_prop2, est_effect_prop)
          weight2 <- c(weight2, weight)
          sim2 <- c(sim2, sim)
          pred_id2 <- c(pred_id2, pred_id)
        }
      }
      
    } else {
      #For models without pO2' term, conditional effect will be 0--add 100 iterations of this for each prediction point
      est_effect_prop2 <- rep(0, 100)
      weight2 <- rep(weights[g], 100)
      pred_id2 <- rep(pred_id, 100)
      sim2 <- rep(1:100, each=1)
    }
    if(g==1){
    df <- data.frame(est_effect_prop=est_effect_prop2, weight=weight2, sim=sim2, pred_id=pred_id2)
    }
    if(g>1){
    df <- bind_rows(df, data.frame(est_effect_prop=est_effect_prop2, weight=weight2, sim=sim2, pred_id=pred_id2))
    }
  }     
  
  #Calculate weighted average
  ens_preds <- df %>% 
    group_by(pred_id, sim) %>% 
    summarise(weighted_mean=weighted.mean(est_effect_prop,weight, na.rm=T))  %>% 
    ungroup()
  
  #Ensemble w/ SD
ens_preds2 <- ens_preds %>% 
    group_by(pred_id) %>% 
    summarise(ensemble_mean=mean(weighted_mean, na.rm=T), ensemble_sd=sd(weighted_mean, na.rm=T)) %>%
    ungroup()
  ens_preds2 <- ens_preds2 %>% 
    mutate(ensemble_lower=(ensemble_mean-ensemble_sd), ensemble_upper=(ensemble_mean+ensemble_sd))%>% 
    ungroup()
  return(ens_preds2)
}

process_species <- function(species2do){
  this_species <- species2do
  bp_est3 <- filter(bp_est2, species==this_species)
  this_data <- bp_est3$data
 
  print(this_species)
  #Calculate metabolic index from correct taxa parameters
  #Pull correct 
  species_tab <- filter(species_table, common_name==this_species)
  #Calculate metabolic index from correct taxa parameters
  this_taxa <- species_tab$MI_Taxa[species_tab$common_name==this_species]
  mi_pars2 <- filter(mi_pars,Group==this_taxa)
  
  Eo <- c(mi_pars2$Eolow, mi_pars2$Eo, mi_pars2$Eohigh)
  
  #New dataframe
  preds <- grid
  
  #Filter by depth and latitude
  #Depth buffer (200m than typical depth, and 0.5 degree latitude north and south)
  #Filter depth and latitude for species
  preds <- filter(preds, depth_m<(species_tab$depth+200))
  southern_limit <- species_tab$southern_limit
  if(!is.na(southern_limit)){
    preds <- filter(preds, latitude>southern_limit)
  }
  
  #Make pred id
  preds$pred_id <- 1:nrow(preds)
  
  #Calculate metabolic index for species
  preds$mi1 <- calc_mi(po2=preds$po2, inv.temp=preds$invtemp, Eo=Eo[1],fancy=F)
  preds$mi2 <- calc_mi(po2=preds$po2, inv.temp=preds$invtemp, Eo=Eo[2],fancy=F)
  preds$mi3 <- calc_mi(po2=preds$po2, inv.temp=preds$invtemp, Eo=Eo[3],fancy=F)
  
  #Pull the data file for scaling
  this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_data, "_dat.rds")))
  #Mean and sd for unscaling MIs
  mean_mi1 <- mean(this_datframe$mi1, na.rm=T)
  sd_mi1 <- sd(this_datframe$mi1, na.rm=T)
  mean_mi2 <- mean(this_datframe$mi2, na.rm=T)
  sd_mi2 <- sd(this_datframe$mi2, na.rm=T)
  mean_mi3 <- mean(this_datframe$mi3, na.rm=T)
  sd_mi3 <- sd(this_datframe$mi3, na.rm=T)
  
  #Calculate
  preds$mi1_s <- (preds$mi1-mean_mi1)/sd_mi1
  preds$mi2_s <- (preds$mi2-mean_mi2)/sd_mi2
  preds$mi3_s <- (preds$mi3-mean_mi3)/sd_mi3
  mi1_s <- preds$mi1_s
  mi2_s <- preds$mi2_s
  mi3_s <- preds$mi3_s
  pred_id <- preds$pred_id

  #Which aic/data to pull
  this_aic <- filter(aic, species==this_species& `data type`==this_data)
  this_dat <- this_data
  print(this_dat)
  #Just the first 5 columns
  this_aic <- select(this_aic, model_names)
  #Only the columns that are not NAs (which means they didn't pass sanity checks)
  #Remove columns that are NAs
  this_aic <- this_aic %>%
    select(where(~ !any(is.na(.))))

  #Get list of these columns (these are the model fits to pull for model averaging)
  models <- colnames(this_aic)
  
  #Pull the model output files
  model_fits <- list()
  for(h in 1:length(models)){
    fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", models[h], ".rds")))
    model_fits[[h]] <- fit
  }
  #Get model weights
  aics <- list()
  for (k in 1:length(models)){
    aic_models <- stats::AIC(model_fits[[k]])
    aics[[k]] <- aic_models
  }
  aics <- unlist(aics)
  delta_aic <- aics - min(aics)
  weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
  set.seed(459384) # for reproducibility and consistency
  #Get lists of parameters to use
  pars_list <- list()
  for(g in 1:length(model_fits)){
    if(models[g] %in% mi_models2use){
      fit <- model_fits[[g]]
      # Get variance-covariance matrix
      V <- fit[["sd_report"]][["cov.fixed"]]
      #Pull out just threshold
      vcov <- V[grep("threshold", rownames(V)), grep("threshold", colnames(V))]
      slope <- fit[["sd_report"]][["value"]][["s_slope"]]
      thresh <- fit[["sd_report"]][["value"]][["s_cut"]]
      #Pull 100 par combinations from variance-covariance matrix
      bp_pars <- MASS::mvrnorm(n=100, mu=c(log(slope),thresh), Sigma=vcov)
      slope_pars <- exp(bp_pars[,1])
      thresh_pars <- bp_pars[,2]
      par <- list(thresh_pars=thresh_pars, slope_pars=slope_pars)
      pars_list[[g]] <- par
    } else {
      pars_list[[g]] <- NA
    }
  }

  #Calculate for each grid cell, parallelized
  #use mcmapply to parallelize
  n <- nrow(preds)
  results <- mcmapply(
    FUN = function(z) {
      calc_grid(
        mi1_s = mi1_s[z],
        mi2_s = mi2_s[z],
        mi3_s = mi3_s[z],
        pred_id=pred_id[z],
        pars_list = pars_list,
        models = models,
        model_fits = model_fits,
        mi_models2use = mi_models2use,
        weights = weights,
        n=n
      )
    },
    z = 1:n,
    mc.cores = 8,
    SIMPLIFY = FALSE
  )
  #convert results to single dataframe
  preds5 <- as.data.frame(bind_rows(results))

  #Add prediction dataframe back
  preds5 <- left_join(preds, preds5, by="pred_id")
  #Other columns
  preds5$data <- this_data 
  preds5$common_name <- this_species
  preds5$depth_limit <- species_tab$depth
  preds5$min_lat <- species_tab$southern_limit
  if(!grepl("iphc", this_dat)){
    preds5$data_type <- "bottom trawl only"
  }
  if(grepl("iphc", this_dat)){
    preds5$data_type <- "bottom trawl & IPHC"
  }
  saveRDS(preds5, file = paste0("output/", output_folder, "/", this_species, "_", "grid.rds"))
  return(preds5)
}

#Apply for all the species
use_previous <- T
if(!use_previous){
species_list <- filter(bp_est2, grepl("coastwide|bc|cc", data))
species_list <- unique(species_list$species)
#apply to all species (if need to run)
grids <- list()
for(i in 1:10){
grids_x <- process_species(species_list[i])
grids[[i]] <- grids_x
}
}

#Or load in
if(use_previous){
grids <- list()
species_list <- filter(bp_est2, grepl("coastwide|bc|cc", data))
species_list <- unique(species_list$species)
for(i in 1:length(species_list)){
  this_species <- species_list[i]
  grid_x <- readRDS(paste0("output/", output_folder, "/", this_species, "_", "grid.rds"))
  grids[[i]] <- grid_x
}
}

#Filter to one year, and to the appropriate region
for(i in 1:length(grids)){
grids2use <- grids[[i]]
grids2use <- filter(grids2use, year==2021)
#grids2use <- filter(grids2use,depth_m<(depth_limit+200))
#grids2use <- filter(grids2use, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))
this_dat <- bp_est2$data[bp_est2$species==species_list[i]]
if(grepl("bc", this_dat)){
  grids2use <- filter(grids2use, region=="bc")
}
if(grepl("cc", this_dat)){
  grids2use <- filter(grids2use, region=="cc")
}

#Remove CC points in Bering Skate, because basically not there
if(species_list[i]=="bering skate"){
  grids2use <- filter(grids2use, !(region=="cc"))
}

grids2use$data_type <- this_dat

if(i==1){
  grids2 <- grids2use
}
if(i>1){
  grids2 <- bind_rows(grids2use, grids2)
}
}

theme_set(theme_bw(base_size = 25))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

##Plot one year, each species, full region
#Plot coastwide models
ylims <- c(min(grids2$Y)*1000, max(grids2$Y)*1000)
xlims <- c(min(grids2$X)*1000, max(grids2$X)*1000)

ggplot(us_coast_proj) + geom_sf() +
  geom_point(filter(grids2, ensemble_mean==0&grepl("coastwide", data_type)),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.1, alpha=0.1)+
  geom_point(filter(grids2,ensemble_mean>0&grepl("coastwide", data_type)),mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.1)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-130,-120), limits=c(xlims))+
  ylim(ylims)+
  facet_wrap("common_name", ncol=5, labeller=labeller(common_name=label_wrap_gen(10)))+
  theme_minimal(base_size=16)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top", panel.spacing = unit(1.5, "lines"), strip.text=element_text(size=13))+
  scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_combined_2021_coastwide.png"),
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

#Plot one year, each species, full region
grids2b <- filter(grids2, region=="bc")

ylims <- c(min(grids2b$Y)*1000, max(grids2b$Y)*1000)
xlims <- c(min(grids2b$X)*1000, max(grids2b$X)*1000)

ggplot(us_coast_proj) + geom_sf() +
  geom_point(filter(grids2b, ensemble_mean==0&grepl("bc", data_type)),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.1, alpha=0.1)+
  geom_point(filter(grids2b,ensemble_mean>0&grepl("bc", data_type)),mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.1)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-130,-120), limits=c(xlims))+
  ylim(ylims)+
  facet_wrap("common_name", ncol=5, labeller=labeller(common_name=label_wrap_gen(10)))+
  theme_minimal(base_size=16)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top", panel.spacing = unit(1.5, "lines"), strip.text=element_text(size=13))+
  scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_combined_2021_bc.png"),
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
##Just California Current
grids2a <- filter(grids2, region=="cc")

ylims <- c(min(grids2a$Y)*1000, max(grids2a$Y)*1000)
xlims <- c(min(grids2a$X)*1000, max(grids2a$X)*1000)

ggplot(us_coast_proj) + geom_sf() +
  geom_point(filter(grids2a, ensemble_mean==0&grepl("cc", data_type)),mapping=aes(x=X*1000, y=Y*1000), colour="grey", size=0.1, alpha=0.1)+
  geom_point(filter(grids2a,ensemble_mean>0&grepl("cc", data_type)),mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.1)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-130,-120), limits=c(xlims))+
  ylim(ylims)+
  facet_wrap("common_name", ncol=5, labeller=labeller(common_name=label_wrap_gen(10)))+
  theme_minimal(base_size=16)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top", panel.spacing = unit(1.5, "lines"), strip.text=element_text(size=13))+
  scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_combined_2021_cc.png"),
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

#Plot one year, each species, WA state
###WASHINGTON COAST
#Plot just one year of data, all species on same plot, zoom to WA
dat2plot_wa <- filter(grids2, latitude>46&latitude<48.5)
#Just species with a CC or coastwide model
dat2plot_wa <- filter(dat2plot_wa, grepl("cc", data_type)|grepl("coastwide", data_type))
#Remove Bering skate, because so little of its range is in the California Current
dat2plot_wa <- filter(dat2plot_wa, !(common_name=="bering skate"))

ggplot(us_coast_proj) + geom_sf() +
  geom_point(dat2plot_wa,mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-125,-120), limits=c(min(dat2plot_wa$X)*1000, max(dat2plot_wa$X)*1000))+
  ylim(min(dat2plot_wa$Y)*1000, max(dat2plot_wa$Y)*1000)+
  facet_wrap("common_name", ncol=5, labeller=labeller(common_name=label_wrap_gen(10)))+
  theme_minimal(base_size=18)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position=c(0.94,0.1), panel.spacing=unit(1,"lines"))+
  scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_combined_2021_WA.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 9,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Compare two years, in WA state
for(i in 1:length(grids)){
  grids2use <- grids[[i]]
  grids2use <- filter(grids2use, year==2010|year==2021)
#  grids2use <- filter(grids2use,depth_m<(depth_limit+200))
#  grids2use <- filter(grids2use, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))
  grids2use <- filter(grids2use, latitude>46&latitude<48.5)
  this_dat <- bp_est2$data[bp_est2$species==species_list[i]]
  grids2use$data_type <- this_dat
  if(i==1){
    grids_b <- grids2use
  }
  if(i>1){
    grids_b <- bind_rows(grids2use, grids_b)
  }
}

#Plot all years, each species separately
for(i in 1:length(grids)) {
  to_do <- grids[[i]]
  this_species <- unique(to_do$common_name)
  print(this_species)
  this_dat <- unique(to_do$data)
 # to_do <- filter(to_do, depth_m<(depth_limit+200))
# to_do<- filter(to_do, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(to_do,mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    scale_x_continuous(breaks=c(-125,-120), limits=c(min(to_do$X)*1000, max(to_do$X)*1000))+
    ylim(min(to_do$Y)*1000, max(to_do$Y)*1000)+
    facet_wrap("year", ncol=5)+
    theme_minimal(base_size=16)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(legend.position=c(0.8,0.1))+
    #ggtitle(paste(this_species, this_dat, sep=" "))
    scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

  ggsave(
    paste0("output/", output_folder, "/", "plots/grid_annual/map_",this_species,"_", this_dat,".png"),
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

#Plot all years, each species separately, WA only
for(i in 1:length(grids)) {
  to_do <- grids[[i]]
  this_species <- unique(to_do$common_name)
  print(this_species)
  this_dat <- unique(to_do$data)
  #to_do <- filter(to_do, depth_m<(depth_limit+200))
  #to_do<- filter(to_do, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))
  to_do <- filter(to_do, latitude>46&latitude<48.5)
  
  ggplot(us_coast_proj) + geom_sf() +
    geom_point(to_do,mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
    #geom_point(filter(dat2plot, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
    scale_x_continuous(breaks=c(-125,-120), limits=c(min(to_do$X)*1000, max(to_do$X)*1000))+
    ylim(min(to_do$Y)*1000, max(to_do$Y)*1000)+
    facet_wrap("year", ncol=5)+
    theme_minimal(base_size=16)+
    xlab("Longitude")+
    ylab("Latitude")+
    theme(legend.position=c(0.8,0.1))+
    #ggtitle(paste(this_species, this_dat, sep=" "))+
    scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")
  
  ggsave(
    paste0("output/", output_folder, "/", "plots/grid_annual_WA/map_",this_species,"_", this_dat,".png"),
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

#Annual proportion of habitat by year
for(i in 1:length(grids)){
test <- grids[[i]]
#test <- filter(test, depth_m<(depth_limit+200))
#test <- filter(test, (latitude>(min_lat-0.5))&(latitude<max_lat+0.5))
#Filter out BC or CC if not that region
this_dat <- bp_est2$data[bp_est2$species==species_list[i]]
if(grepl("bc", this_dat)){
  test <- filter(test, region=="bc")
}
if(grepl("cc", this_dat)){
  test <- filter(test, region=="cc")
}

#Remove CC points in Bering Skate, because basically not there
if(species_list[i]=="bering skate"){
  test<- filter(test, !(region=="cc"))
}

area_sum2 <- filter(test, year==2021)
areas <- area_sum2 %>%
  group_by(region) %>%
  summarise(area_sum=sum(area))  %>%
  ungroup()

test <- left_join(test, areas, by=c("region"))

test2 <- filter(test, ensemble_mean>0.1)
test3 <- filter(test, ensemble_lower>0.1)
test4 <- filter(test, ensemble_upper>0.1)

effects_full_annual_summary <- test2 %>%
  group_by(year, region, common_name) %>%
  summarise(area_sum=sum(area)/area_sum) %>%
  distinct()%>%
  ungroup()

effects_full_annual_summary2 <- test3 %>%
  group_by(year, region, common_name) %>%
  summarise(area_sum=sum(area)/area_sum) %>%
  distinct()%>%
  ungroup()
effects_full_annual_summary3 <- test4 %>%
  group_by(year, region, common_name) %>%
  summarise(area_sum=sum(area)/area_sum) %>%
  distinct()%>%
  ungroup()

effects_full_summary <- full_join(effects_full_annual_summary, effects_full_annual_summary2, by=c("year", "region", "common_name"))
effects_full_summary$se1 <- effects_full_summary$area_sum.y
effects_full_summary$area_sum_med <- effects_full_summary$area_sum.x
effects_full_summary$area_sum.y <- NULL
effects_full_summary$area_sum.x <- NULL
effects_full_summary <- full_join(effects_full_summary, effects_full_annual_summary3, by=c("year", "region", "common_name"))
effects_full_summary$se2 <- effects_full_summary$area_sum
effects_full_summary$area_sum <- NULL
effects_full_summary$species <- unique(test$common_name)

if(i==1){
effects <- effects_full_summary
}
if(i>1){
effects <- bind_rows(effects, effects_full_summary)
}
}

#Regions
effects$region <- factor(effects$region, levels=c("bc", "cc"))
labs <- c("British Columbia", "California Current")
names(labs) <- c("bc", "cc")

#plot
ggplot(effects, aes(x=year, y=area_sum_med))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", scales="free_y",ncol=4,labeller=labeller(common_name=label_wrap_gen(15)))+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank(), text=element_text(size=20),
  legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_fill_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  scale_colour_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  xlab("Year") +
  ylab("Proportion of area w/ >10% local density reduction")

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_area_below_grid.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 12,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Just selected species, for subset figure
temp <- filter(grids_b, common_name=="pacific cod"|common_name=="shortspine thornyhead"|common_name=="spotted ratfish"|common_name=="silvergray rockfish")

a <- ggplot(us_coast_proj) + geom_sf() +
  geom_point(temp,mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean), size=0.5)+
  #geom_point(filter(dat2plot_wa, est_effect_raw==max_effect),mapping=aes(x=X*1000, y=Y*1000), colour="#440154FF", size=0.5)+
  #geom_point(filter(dat2plot_wa, est_effect_raw<max_effect),mapping=aes(x=X*1000, y=Y*1000, colour=est_effect_prop), size=0.5)+
  scale_x_continuous(breaks=c(-125,-120), limits=c(min(temp$X)*1000, max(temp$X)*1000))+
  ylim(min(temp$Y)*1000, max(temp$Y)*1000)+
  facet_nested(~common_name~year, labeller=labeller(common_name=label_wrap_gen(15)))+
  theme_minimal(base_size=16)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top", legend.text = element_text(angle = 45, hjust = 1), legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

temp <- filter(effects, common_name=="pacific cod"|common_name=="shortspine thornyhead"|common_name=="spotted ratfish"|common_name=="silvergray rockfish")

b <- ggplot(temp, aes(x=year, y=area_sum_med))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", ncol=1,scales="free_y",labeller=labeller(common_name=label_wrap_gen(15)))+
  theme(legend.title=element_blank(), strip.background = element_blank(), strip.text=element_blank(),
        legend.justification="center", legend.box.spacing = unit(0, "pt"), legend.position="top",
        legend.text=element_text(size=16), axis.text=(element_text(size=16)), axis.title=element_text(size=16))+
  scale_fill_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  scale_colour_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  guides(fill = guide_legend(nrow = 2))+
  xlab("Year") +
  ylab("Proportion of area w/ >10% local density reduction")

patchwork <- b+a
patchwork + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18))

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_area_below_grid_subset.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  width = 8.5,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#####Annual proportion of habitat by depth
for(i in 1:length(grids)){
  test <- grids[[i]]
 #test <- filter(test, depth_m<(depth_limit+200))
  #test <- filter(test, (latitude>(min_lat-0.5))&(latitude<max_lat+0.5))
  #Filter out BC or CC if not that region
  this_dat <- bp_est2$data[bp_est2$species==species_list[i]]
  if(grepl("bc", this_dat)){
    test <- filter(test, region=="bc")
  }
  if(grepl("cc", this_dat)){
    test <- filter(test, region=="cc")
  }
  
  #Remove CC points in Bering Skate, because basically not there
  if(species_list[i]=="bering skate"){
    test<- filter(test, !(region=="cc"))
  }
  
test$depth_bin <- cut(test$depth_m, labels=FALSE, breaks=seq(0, max(test$depth_m), by=25), include.lowest = TRUE)
test$depth_bin <- test$depth_bin*25

test <- filter(test, year==2021)

areas <- test %>%
  group_by(depth_bin, region) %>%
  summarise(total_area=sum(area))  %>%
ungroup()

test <- left_join(test, areas, by=c("region", "depth_bin"))

test2 <- filter(test, ensemble_mean>0.1)
test3 <- filter(test, ensemble_lower>0.1)
test4 <- filter(test, ensemble_upper>0.1)

effects_full_annual_summary <- test2 %>%
  group_by(depth_bin, region, total_area) %>%
  summarize(area_sum = sum(area)) %>%
  mutate(prop=area_sum/total_area) %>%
  ungroup()
effects_full_annual_summary2 <- test3 %>%
  group_by(depth_bin, region, total_area) %>%
  summarise(area_sum = sum(area)) %>%
  mutate(prop=area_sum/total_area) %>%
  ungroup()
effects_full_annual_summary3 <- test4 %>%
  group_by(depth_bin, region, total_area) %>%
  summarise(area_sum = sum(area)) %>%
  mutate(prop=area_sum/total_area) %>%
  ungroup()

effects_full_summary <- full_join(effects_full_annual_summary, effects_full_annual_summary2, by=c("depth_bin", "region"))
effects_full_summary$se1 <- effects_full_summary$prop.y
effects_full_summary$area_sum_med <- effects_full_summary$prop.x
effects_full_summary$prop.y <- NULL
effects_full_summary$prop.x <- NULL
effects_full_summary <- full_join(effects_full_summary, effects_full_annual_summary3, by=c("depth_bin", "region"))
effects_full_summary$se2 <- effects_full_summary$prop
effects_full_summary$prop <- NULL
effects_full_summary$species <- unique(test$common_name)
effects_full_summary$depth_limit <- unique(test$depth_limit)


if(i==1){
  effects2 <- effects_full_summary
}
if(i>1){
  effects2 <- bind_rows(effects2, effects_full_summary)
}
}

#plot
ggplot(effects2, aes(x=depth_bin, y=area_sum_med))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  #geom_vline(mapping=(aes(data=filter(species_table, common_name %in% species_list), xintercept=depth), linetype="dashed")+
  facet_wrap("species",scales="free_x", ncol=5, labeller=labeller(species=label_wrap_gen(15)))+
  theme(legend.position="top")+
ggh4x::facetted_pos_scales(x=list(scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), scale_x_continuous(limits=c(0,800), breaks=c(100,300,500, 700)), NULL, scale_x_continuous(breaks=c(100,300,500, 700), limits=c(0,600)), 
                             scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)),scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), scale_x_continuous(breaks=c(100,300,500, 700), limits=c(0,700)), 
                             scale_x_continuous(breaks=c(100,300,500), limits=c(0,600)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,600)), 
                             scale_x_continuous(breaks=c(100,300,500, 700), limits=c(0,700)), scale_x_continuous(breaks=c(100,500,900), limits=c(0,1200)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,500)), 
                             scale_x_continuous(breaks=c(100,300,500), limits=c(0,600)), scale_x_continuous(breaks=c(100,300,500), limits=c(0,600))))+
  theme(legend.title=element_blank(), strip.background = element_blank(), text=element_text(size=16),
  legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_fill_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  scale_colour_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  xlab("Depth (m)") +
  ylab("Proportion of area w/ >10% local density reduction")

ggsave(
  paste0("output/", output_folder, "/", "plots/grid_area_depth_bins.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8,5,
  height = 10,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Annual with states
#Add states
for(i in 1:length(grids)){
  test <- grids[[i]]
#  test <- filter(test, depth_m<(depth_limit+200))
 # test <- filter(test, (latitude>(min_lat-0.5))&(latitude<max_lat+0.5))
  #add state
  #Filter out BC or CC if not that region
  this_dat <- bp_est2$data[bp_est2$species==species_list[i]]
  if(grepl("bc", this_dat)){
    test <- filter(test, region=="bc")
  }
  if(grepl("cc", this_dat)){
    test <- filter(test, region=="cc")
  }
  
  #Remove CC points in Bering Skate, because basically not there
  if(species_list[i]=="bering skate"){
    test<- filter(test, !(region=="cc"))
  }
  test$state <- case_when(test$latitude<42~"California",
                          test$latitude >42 & test$latitude<46~"Oregon",
                          test$latitude>46& test$region=="cc"~"Washington",
                          test$latitude>46& test$region=="bc"~"British Columbia")
  
  
  area_sum2 <- filter(test, year==2021)
  areas <- area_sum2 %>%
    group_by(state) %>%
    summarise(total_area=sum(area))  %>%
    ungroup()
  
  test <- left_join(test, areas, by=c("state"))
  
  test2 <- filter(test, ensemble_mean>0.1)
  test3 <- filter(test, ensemble_lower>0.1)
  test4 <- filter(test, ensemble_upper>0.1)
  
  effects_full_annual_summary <- test2 %>%
    group_by(year, state, common_name) %>%
    summarise(area_sum=sum(area)/total_area) %>%
    distinct()%>%
    ungroup()
  
  effects_full_annual_summary2 <- test3 %>%
    group_by(year, state, common_name) %>%
    summarise(area_sum=sum(area)/total_area) %>%
    distinct()%>%
    ungroup()
  effects_full_annual_summary3 <- test4 %>%
    group_by(year, state, common_name) %>%
    summarise(area_sum=sum(area)/total_area) %>%
    distinct()%>%
    ungroup()
  
  effects_full_summary <- full_join(effects_full_annual_summary, effects_full_annual_summary2, by=c("year", "state", "common_name"))
  effects_full_summary$se1 <- effects_full_summary$area_sum.y
  effects_full_summary$area_sum_med <- effects_full_summary$area_sum.x
  effects_full_summary$area_sum.y <- NULL
  effects_full_summary$area_sum.x <- NULL
  effects_full_summary <- full_join(effects_full_summary, effects_full_annual_summary3, by=c("year", "state", "common_name"))
  effects_full_summary$se2 <- effects_full_summary$area_sum
  effects_full_summary$area_sum <- NULL
  effects_full_summary$species <- unique(test$common_name)
  
  if(i==1){
    effects3 <- effects_full_summary
  }
  if(i>1){
    effects3 <- bind_rows(effects3, effects_full_summary)
  }
}

#plot
effects3$state <- factor(effects3$state, levels=c("British Columbia", "Washington", "Oregon", "California"))

ggplot(effects3, aes(x=year, y=area_sum_med))+
  geom_line(aes(colour=state))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=state), alpha=0.2)+
  facet_wrap("common_name", scales="free_y",ncol=5,labeller=labeller(common_name=label_wrap_gen(18)))+
  theme(legend.position="top")+
  theme(legend.title=element_blank(), strip.background = element_blank())+
  theme(text=element_text(size=16),   legend.justification="center",
        legend.box.spacing = unit(0, "pt"))+
  scale_fill_manual(values=c("red", "darkgreen", "royalblue","gold"), drop=FALSE)+
  scale_colour_manual(values=c("red", "darkgreen", "royalblue","gold"), drop=FALSE)+
  scale_x_continuous(breaks=c(2008,2016,2024))+
  xlab("Year") +
  ylab("Proportion of area w/ >10% local density reduction")

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_area_below_grid_states.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width =10,
  height = 9,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

###Add full density predictions and density predictions with O2 above threshold
species_list <- filter(bp_est2, grepl("coastwide|bc|cc", data))
species_list <- unique(species_list$species)
grids2 <- list()
#Calculate predictions for each grid cell
for(i in 1:length(species_list)){
  this_species <- species_list[i]
  bp_est3 <- filter(bp_est2, species==this_species)
  this_data <- bp_est3$data
  print(this_species)
  #Calculate metabolic index from correct taxa parameters
  #Pull correct 
  species_tab <- filter(species_table, common_name==this_species)
  #Calculate metabolic index from correct taxa parameters
  this_taxa <- species_tab$MI_Taxa[species_tab$common_name==this_species]
  mi_pars2 <- filter(mi_pars,Group==this_taxa)
  
  Eo <- c(mi_pars2$Eolow, mi_pars2$Eo, mi_pars2$Eohigh)
  
  #New dataframe
  preds <- grids[[i]]
  
  #Filter by depth and latitude
  #Depth buffer (200m than typical depth, and 0.5 degree latitude north and south)
  #Filter depth and latitude for species
  preds <- filter(preds, depth_m<(species_tab$depth+200))
  southern_limit <- species_tab$southern_limit
  if(!is.na(southern_limit)){
    preds <- filter(preds, latitude>southern_limit)
  }
  
  #If BC-only model, filter to BC region
  if(grepl("bc", this_data)){
    preds <- filter(preds, region=="bc")
  }
  
  #If CC-only model, filter to CC regio
  if(grepl("cc", this_data)){
    preds <- filter(preds, region=="cc")
  }
  
  if(this_species=="bering skate"){
    preds <- filter(preds, region=="bc")
  }
  
  #Calculate metabolic index for species
  preds$mi1 <- calc_mi(po2=preds$po2, inv.temp=preds$invtemp, Eo=Eo[1],fancy=F)
  preds$mi2 <- calc_mi(po2=preds$po2, inv.temp=preds$invtemp, Eo=Eo[2],fancy=F)
  preds$mi3 <- calc_mi(po2=preds$po2, inv.temp=preds$invtemp, Eo=Eo[3],fancy=F)
  
  #Pull the data file for scaling
  this_datframe <- try(readRDS(file = paste0("output/", output_folder, "/", this_species, "_", this_data, "_dat.rds")))
  #Mean and sd for unscaling MIs
  mean_mi1 <- mean(this_datframe$mi1, na.rm=T)
  sd_mi1 <- sd(this_datframe$mi1, na.rm=T)
  mean_mi2 <- mean(this_datframe$mi2, na.rm=T)
  sd_mi2 <- sd(this_datframe$mi2, na.rm=T)
  mean_mi3 <- mean(this_datframe$mi3, na.rm=T)
  sd_mi3 <- sd(this_datframe$mi3, na.rm=T)
  
  #Calculate
  preds$mi1_s <- (preds$mi1-mean_mi1)/sd_mi1
  preds$mi2_s <- (preds$mi2-mean_mi2)/sd_mi2
  preds$mi3_s <- (preds$mi3-mean_mi3)/sd_mi3
  mi1_s <- preds$mi1_s
  mi2_s <- preds$mi2_s
  mi3_s <- preds$mi3_s
  
  #Scale depth and temperature
  preds$log_depth <- preds$depth_ln
  preds$log_depth_scaled <- (preds$log_depth - mean(log(this_datframe$depth), na.rm=T))/sd(log(this_datframe$depth), na.rm=T)
  preds$log_depth_scaled2 <- preds$log_depth_scaled^2
  preds$log_depth_scaled3 <- preds$log_depth_scaled^3
  preds$temp_scaled <- (preds$temperature_C - mean(this_datframe$temperature_C, na.rm=T))/sd(this_datframe$temperature_C, na.rm=T)
  preds$temp_scaled2 <- preds$temp_scaled^2
  
  #add correct survey name
  preds$survey <- ifelse(preds$region=="cc", "nwfsc", "dfo")
  
  ##Filter to only years of data
  years2use <- unique(this_datframe$year)
  
  #Filter
  preds <- filter(preds, year %in% years2use)
  preds$year <- factor(preds$year)
  
  #Prediction ID (row number)
  preds$pred_id <- 1:nrow(preds)
  
  #Which aic/data to pull
  this_aic <- filter(aic, species==this_species& `data type`==this_data)
  this_dat <- this_data
  print(this_dat)
  #Just the first 5 columns
  this_aic <- select(this_aic, model_names)
  #Only the columns that are not NAs (which means they didn't pass sanity checks)
  #Remove columns that are NAs
  this_aic <- this_aic %>%
    select(where(~ !any(is.na(.))))
  
  #Get list of these columns (these are the model fits to pull for model averaging)
  models <- colnames(this_aic)
  
  #Pull the model output files
  model_fits <- list()
  for(h in 1:length(models)){
    fit <- try(readRDS(file = paste0("output/", output_folder,"/", this_species,"_",this_dat,"_", models[h], ".rds")))
    model_fits[[h]] <- fit
  }
  #Get model weights
  aics <- list()
  for (k in 1:length(models)){
    aic_models <- stats::AIC(model_fits[[k]])
    aics[[k]] <- aic_models
  }
  aics <- unlist(aics)
  delta_aic <- aics - min(aics)
  weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
  set.seed(459384) # for reproducibility and consistency
  
  for(g in 1:length(model_fits)){
    #Calculate the breakpoint effect
    #Prediction
    print(models[[g]])
    prediction <- predict(fit, newdata=preds, se.fit=FALSE)
    prediction$weight <- weights[g]
    if(g==1){
      df <- prediction
    }
    if(g>1){
      df <- bind_rows(df, prediction)
    }
  }     
  
  #Calculate weighted average
  ens_preds <- df %>% 
    group_by(pred_id) %>% 
    summarise(weighted_mean=weighted.mean(est,weight, na.rm=T))  %>% 
    ungroup()
  
  #Ensemble w/ SD
  ens_preds2 <- ens_preds %>% 
    group_by(pred_id) %>% 
    summarise(ensemble_mean=mean(weighted_mean, na.rm=T))
  
  #Add prediction dataframe back
  preds5 <- left_join(preds, ens_preds2, by="pred_id")
  
  #And for super high O2
  #Add dummy MI columns for MI above threshold
  #Keep real values
  preds$mi1_s_real <- preds$mi1_s
  preds$mi2_s_real <- preds$mi2_s
  preds$mi3_s_real <- preds$mi3_s
  preds$mi1_s <- 1000
  preds$mi2_s <- 1000
  preds$mi3_s <- 1000
  
  for(g in 1:length(model_fits)){
    #Calculate the breakpoint effect
    #Prediction
    prediction <- predict(fit, newdata=preds, se.fit=FALSE)
    prediction$weight <- weights[g]
    if(g==1){
      df <- prediction
    }
    if(g>1){
      df <- bind_rows(df, prediction)
    }
  }     
  
  #Calculate weighted average
  ens_preds <- df %>% 
    group_by(pred_id) %>% 
    summarise(weighted_mean=weighted.mean(est,weight, na.rm=T))  %>% 
    ungroup()
  
  #Ensemble w/ SD
  ens_preds2 <- ens_preds %>% 
    group_by(pred_id) %>% 
    summarise(ensemble_mean_full=mean(weighted_mean, na.rm=T))
  
  #Add back
  #Add prediction dataframe back
  preds5 <- left_join(preds5, ens_preds2, by="pred_id")
  
  #Other columns
  preds5$data <- this_data 
  preds5$common_name <- this_species
  preds5$depth_limit <- species_tab$depth
  preds5$min_lat <- species_tab$southern_limit
  
  #Make MI back to normal
  preds5$mi1_s <- preds$mi1_s_real
  preds5$mi2_s <- preds$mi2_s_real
  preds5$mi3_s <- preds$mi3_s_real
  #remove fake MI
  preds5$mi1_s_real <- NULL
  preds5$mi2_s_real <- NULL
  preds5$mi3_s_real <- NULL
  
  if(!grepl("iphc", this_dat)){
    preds5$data_type <- "bottom trawl only"
  }
  if(grepl("iphc", this_dat)){
    preds5$data_type <- "bottom trawl & IPHC"
  }
  preds5$dat <- this_data
  saveRDS(preds5, file = paste0("output/", output_folder, "/", this_species, "_", "grid2.rds"))
  grids2[[i]] <- preds5
}

#Add species names to list to label each dataset
names(grids2) <- species_list

##Just 2010 and 2021 for plotting
for(i in 1:length(grids2)){
  grids2use <- grids2[[i]]
  grids2use <- filter(grids2use, year==2010|year==2021)
  #  grids2use <- filter(grids2use,depth_m<(depth_limit+200))
  #  grids2use <- filter(grids2use, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))
  grids2use <- filter(grids2use, latitude>46&latitude<48.5)
  #Scale predicted biomass
  grids2use$density_scaled <- exp(grids2use$ensemble_mean_full)/exp(max(grids2use$ensemble_mean_full))
  if(i==1){
    grids_c <- grids2use
  }
  if(i>1){
    grids_c <- bind_rows(grids2use, grids_c)
  }
}

#Filter to species of interest
temp <- filter(grids_c, common_name=="pacific halibut"|common_name=="shortspine thornyhead"|common_name=="spotted ratfish"|common_name=="pacific cod")

##Plot reduction from oxygen, faded by biomass
a <- ggplot(us_coast_proj) + geom_sf() +
  geom_point(filter(temp,density_scaled>0.05),mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean.x), size=0.1)+
  geom_point(filter(temp,density_scaled<0.05),mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean.x), size=0.1, alpha=0.04)+
  scale_x_continuous(breaks=c(-125,-120), limits=c(min(temp$X)*1000, max(temp$X)*1000))+
  ylim(min(temp$Y)*1000, max(temp$Y)*1000)+
  facet_nested(~common_name~year, labeller=labeller(common_name=label_wrap_gen(15)))+
  theme_minimal(base_size=16)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top", legend.text = element_text(angle = 45, hjust = 1), legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

temp2 <- filter(effects, common_name=="pacific halibut"|common_name=="shortspine thornyhead"|common_name=="spotted ratfish"|common_name=="pacific cod")

b <- ggplot(temp2, aes(x=year, y=area_sum_med))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", ncol=1,scales="free_y",labeller=labeller(common_name=label_wrap_gen(15)))+
  theme(legend.title=element_blank(), strip.background = element_blank(), strip.text=element_blank(),
        legend.justification="center", legend.box.spacing = unit(0, "pt"), legend.position="top",
        legend.text=element_text(size=16), axis.text=(element_text(size=16)), axis.title=element_text(size=16))+
  scale_fill_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  scale_colour_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  guides(fill = guide_legend(nrow = 2))+
  xlab("Year") +
  ylab("Proportion of area w/ >10% local density reduction")

patchwork <- b+a
patchwork + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18))

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_area_below_grid_subset_alpha_biomass.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  width = 8.5,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Plot side-by-side panels of reduction in local density and reductionin biomass
a <- ggplot(us_coast_proj) + geom_sf() +
  geom_point(temp,mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean.x), size=0.1)+
  scale_x_continuous(breaks=c(-125,-120), limits=c(min(temp$X)*1000, max(temp$X)*1000))+
  ylim(min(temp$Y)*1000, max(temp$Y)*1000)+
  facet_nested(~common_name~year)+
  theme_minimal(base_size=16)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(strip.text.y = element_blank(),legend.position="top", legend.text = element_text(angle = 45, hjust = 1),
    legend.direction = "horizontal",
    plot.margin = unit(c(0, 0, 0, 0), "null"),
    legend.title.align = 0.5,
    legend.title.position = "top",
    legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

b <-ggplot(us_coast_proj) + geom_sf() +
  geom_point(temp,mapping=aes(x=X*1000, y=Y*1000, colour=density_scaled), size=0.1)+
  scale_x_continuous(breaks=c(-125,-120), limits=c(min(temp$X)*1000, max(temp$X)*1000))+
  ylim(min(temp$Y)*1000, max(temp$Y)*1000)+
  facet_nested(~common_name~year, labeller=labeller(common_name=label_wrap_gen(15)))+
  theme_minimal(base_size=16)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        legend.position="top", legend.text = element_text(angle = 45, hjust = 1),
              legend.direction = "horizontal",
              legend.title.align = 0.5,
              legend.title.position = "top",
              legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_colour_viridis(name="Total Predicted \n Density (Scaled) ", limits=c(0,1), breaks=c(0,0.25, 0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1), oob = scales::squish, option="plasma")

c <- ggplot(temp2, aes(x=year, y=area_sum_med))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", ncol=1,scales="free_y",labeller=labeller(common_name=label_wrap_gen(15)))+
  theme(legend.title=element_blank(), strip.background = element_blank(), strip.text=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        legend.justification="center", legend.box.spacing = unit(0, "pt"), legend.position="top",
        legend.text=element_text(size=16), axis.text=(element_text(size=16)), axis.title=element_text(size=16))+
  scale_fill_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  scale_colour_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  guides(fill = guide_legend(nrow = 2))+
  xlab("Year") +
  ylab("Proportion of area w/ >10% local density reduction")

p <- c |a |b
#Change proportions
p + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18))


ggsave(
  paste0("output/", output_folder, "/", "plots/prop_area_below_grid_subset_biomass_also.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  width = 16,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Plot with circles of predicted density
#Thin out points for plotting biomass
#Every 0.25 of a degree latitude
temp3 <- temp %>%
  group_by(common_name, year) %>%
  slice(seq(1, n(), by=300)) %>%
  ungroup()

a <- ggplot(us_coast_proj) + geom_sf() +
  geom_point(temp,mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean.x), size=0.1)+
  geom_point(temp3,mapping=aes(x=X*1000, y=Y*1000, size=density_scaled), shape=1)+
  scale_x_continuous(breaks=c(-125,-120), limits=c(min(temp$X)*1000, max(temp$X)*1000))+
  ylim(min(temp$Y)*1000, max(temp$Y)*1000)+
  facet_nested(~common_name~year, labeller=labeller(common_name=label_wrap_gen(15)))+
  theme_minimal(base_size=16)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position="top", 
        legend.text = element_text(angle = 45, hjust = 1), 
        legend.justification="center", 
        legend.box.spacing = unit(0, "pt"),
        legend.direction = "horizontal",
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        legend.title.align = 0.5,
        legend.title.position = "top")+
  scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

temp2 <- filter(effects, common_name=="pacific halibut"|common_name=="shortspine thornyhead"|common_name=="spotted ratfish"|common_name=="pacific cod")

b <- ggplot(temp2, aes(x=year, y=area_sum_med))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", ncol=1,scales="free_y",labeller=labeller(common_name=label_wrap_gen(15)))+
  theme(legend.title=element_blank(), strip.background = element_blank(), strip.text=element_blank(),
        legend.direction = "horizontal",
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        legend.title.align = 0.5,
        legend.title.position = "top",
        legend.justification="center",
        legend.box.spacing = unit(0, "pt"),
        legend.text=element_text(size=16), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=16))+
  scale_fill_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  scale_colour_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  guides(fill = guide_legend(nrow = 2))+
  xlab("Year") +
  ylab("Proportion of area w/ >10% local density reduction")

patchwork <- b+a
patchwork + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18))

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_area_below_grid_subset_circles.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  width = 8.5,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Scatterplot of reduction vs predicted biomass
a <- ggplot(temp, aes(x=density_scaled, y=ensemble_mean.x))+
  geom_point(aes(color=as.factor(year)),alpha=0.5, size=0.2)+
  facet_wrap("common_name", scales="fixed", ncol=1)+
  xlab("Scaled Predicted Density")+
  ylab("Reduction in Local Density")+
  scale_colour_manual(values=c("2010"="#0B040Fff", "2021"="#357BA2FF"), name="Year")+
  theme_minimal(base_size=16)+
  theme(legend.position="top",
legend.direction = "horizontal",
legend.title.align = 0.5,
legend.title.position = "top",
legend.justification="center", legend.box.spacing = unit(0, "pt"),
strip.background = element_blank(), strip.text=element_blank(),)+
  guides(color = guide_legend(override.aes = list(size = 10, alpha=1)))

b <-ggplot(us_coast_proj) + geom_sf() +
  geom_point(temp,mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean.x), size=0.1)+
  scale_x_continuous(breaks=c(-125,-120), limits=c(min(temp$X)*1000, max(temp$X)*1000))+
  ylim(min(temp$Y)*1000, max(temp$Y)*1000)+
  facet_nested(~common_name~year, labeller=labeller(common_name=label_wrap_gen(15)))+
  theme_minimal(base_size=16)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        legend.position="top", legend.text = element_text(angle = 45, hjust = 1),
        legend.direction = "horizontal",
        legend.title.align = 0.5,
        legend.title.position = "top",
        legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

c <- ggplot(temp2, aes(x=year, y=area_sum_med))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", ncol=1,scales="free_y",labeller=labeller(common_name=label_wrap_gen(15)))+
  theme(legend.title=element_blank(), strip.background = element_blank(), strip.text=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        legend.justification="center", legend.box.spacing = unit(0, "pt"), legend.position="top",
        legend.text=element_text(size=16), axis.text=(element_text(size=16)), axis.title=element_text(size=16))+
  scale_fill_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  scale_colour_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  guides(fill = guide_legend(nrow = 2))+
  xlab("Year") +
  ylab("Proportion of area w/ >10% local density reduction")

p <- c |a |b
#Change proportions
p + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18))

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_area_below_grid_subset_scatterplot.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  width = 10,
  height = 11,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Plot the scatterplot of predicted biomass to percentage reduction from O2 for all species
for(i in 1:length(grids2)){
  grids2use <- grids2[[i]]
  grids2use <- filter(grids2use, year==2017)
  #  grids2use <- filter(grids2use,depth_m<(depth_limit+200))
  #  grids2use <- filter(grids2use, (latitude>(min_lat-0.5))&(latitude<(max_lat+0.5)))
  grids2use <- filter(grids2use, latitude>46&latitude<48.5)
  #Scale predicted biomass
  grids2use$density_scaled <- exp(grids2use$ensemble_mean_full)/exp(max(grids2use$ensemble_mean_full))
  if(i==1){
    grids_d <- grids2use
  }
  if(i>1){
    grids_d <- bind_rows(grids2use, grids_d)
  }
}

ggplot(grids_d, aes(x=density_scaled, y=ensemble_mean.x))+
  geom_point(size=0.1, alpha=0.1)+
  facet_wrap("common_name", scales="free_y", ncol=1)+
  xlab("Scaled Predicted Density")+
  ylab("Reduction in Local Density")+
  facet_wrap("common_name", labeller=labeller(common_name=label_wrap_gen(15)))+
 # scale_colour_manual(values=c("2010"="#440154", "2021"="#FDE725FF"), name="Year")+
  theme(legend.position="top",
        legend.direction = "horizontal",
        legend.title.align = 0.5,
        legend.title.position = "top",
        legend.justification="center", legend.box.spacing = unit(0, "pt"))
 # guides(color = guide_legend(override.aes = list(size = 10, alpha=1)))

ggsave(
  paste0("output/", output_folder, "/", "plots/scatterplot_density_reduction_all_species.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  width = 16,
  height = 14,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)


##Plot with predicted density by depth, for a sample latitude (46 deg), for year 2021, with full and actual oxygen
test <- filter(temp, latitude>47&latitude<47.1&year==2021)

#Scale biomass to maximum biomass in that species
test<- test %>%
  group_by(common_name) %>%
  mutate(ensemble_mean.y2= exp(ensemble_mean.y)/max(exp(ensemble_mean_full), na.rm=T),
         ensemble_mean_full2 = exp(ensemble_mean_full)/max(exp(ensemble_mean_full),na.rm=T))  %>%
  ungroup()

##Line plot of predicted vs full biomass across depths
a <-ggplot(test, aes(x=depth_m, y=ensemble_mean.y2))+
  geom_line(aes(colour="Density w/ oxygen reduction"))+
  geom_line(data=test, aes(x=depth_m, y=ensemble_mean_full2,colour="Full Density"))+
  facet_wrap("common_name", ncol=1,scales="free",labeller=labeller(common_name=label_wrap_gen(15)))+
  scale_colour_manual(values=c("#ffcf20FF", "#482577FF"), name="")+
  theme(legend.title=element_blank(), strip.background = element_blank(), strip.text=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        legend.justification="center", legend.box.spacing = unit(0, "pt"), legend.position="top",
        legend.text=element_text(size=25), axis.text=(element_text(size=25)), axis.title=element_text(size=25))+
  guides(colour = guide_legend(nrow = 2))+
  ylab("Scaled Fish Density")+
  xlab("Depth (m)")


b <-ggplot(us_coast_proj) + geom_sf() +
  geom_point(temp,mapping=aes(x=X*1000, y=Y*1000, colour=ensemble_mean.x), size=0.1)+
  scale_x_continuous(breaks=c(-125,-120), limits=c(min(temp$X)*1000, max(temp$X)*1000))+
  scale_y_continuous(limits=c(min(temp$Y)*1000, max(temp$Y)*1000), breaks=c(46,47,48))+
  facet_nested(~common_name~year, labeller=labeller(common_name=label_wrap_gen(15)))+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        legend.position="top", legend.text = element_text(angle = 45, hjust = 1, size=25),
        legend.direction = "horizontal",
        legend.title.align = 0.5,
        legend.title.position="top",
        axis.text=(element_text(size=25)), axis.title=element_text(size=25), 
        strip.text=element_text(size=25),
        legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  scale_colour_viridis(name="Reduction in \nlocal density ", limits=c(0,0.8), breaks=c(0,0.25, 0.5,0.75,1), labels=scales::percent(c(0,0.25,0.5,0.75,1)), oob = scales::squish, option="plasma")

c <- ggplot(temp2, aes(x=year, y=area_sum_med))+
  geom_line(aes(colour=region))+
  geom_ribbon(aes(ymin=se1, ymax=se2, fill=region), alpha=0.2)+
  facet_wrap("common_name", ncol=1,scales="free_y",labeller=labeller(common_name=label_wrap_gen(15)))+
  theme(legend.title=element_blank(), strip.background = element_blank(), strip.text=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        legend.justification="center", legend.box.spacing = unit(0, "pt"), legend.position="top",
        legend.text=element_text(size=25), axis.text=(element_text(size=25)), axis.title=element_text(size=25))+
  scale_fill_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  scale_colour_manual(values=c("#44AA99","#CC6677"), drop=FALSE, labels=labs )+
  guides(fill = guide_legend(nrow = 2))+
  xlab("Year") +
  ylab("Proportion of area w/ >10% local density reduction")

p <- c |a |b

#Change proportions
p + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 18))

ggsave(
  paste0("output/", output_folder, "/", "plots/prop_area_below_grid_subset_line_plot.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  width = 16,
  height = 14,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#All species, line plots
##Plot with predicted density by depth, for a sample latitude (46 deg), for year 2021, with full and actual oxygen
test <- filter(grids_c, latitude>47&latitude<47.1&year==2021)
test <- filter(test, data=="coastwide"|data=="cc")

#Scale biomass to maximum biomass in that species
test<- test %>%
  group_by(common_name) %>%
  mutate(ensemble_mean.y2= exp(ensemble_mean.y)/max(exp(ensemble_mean_full), na.rm=T),
         ensemble_mean_full2 = exp(ensemble_mean_full)/max(exp(ensemble_mean_full),na.rm=T))  %>%
  ungroup()

##Line plot of predicted vs full biomass across depths
ggplot(test, aes(x=depth_m, y=ensemble_mean.y2))+
  geom_line(aes(colour="Density w/ oxygen reduction"))+
  geom_line(data=test, aes(x=depth_m, y=ensemble_mean_full2,colour="Full Density"))+
  facet_wrap("common_name", scales="free",labeller=labeller(common_name=label_wrap_gen(15)))+
  scale_colour_manual(values=c("#ffcf20FF", "#482577FF"), name="")+
  theme(legend.title=element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        legend.justification="center", legend.box.spacing = unit(0, "pt"), legend.position="top",
        legend.text=element_text(size=25), axis.text=(element_text(size=25)), axis.title=element_text(size=25))+
  guides(colour = guide_legend(nrow = 2))+
  ylab("Scaled Fish Density")+
  xlab("Depth (m)")

ggsave(
  paste0("output/", output_folder, "/", "plots/full_reduced_biomass_line_plots.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  width = 16,
  height = 14,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)
