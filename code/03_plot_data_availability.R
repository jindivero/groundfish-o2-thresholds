library(sdmTMB)
library(dplyr)
library(Metrics)
library(ggplot2)
library(tidyr)
library(visreg)
library(ggpubr)
library(purrr)
library(readxl)
library(openxlsx2)
library(ggpattern)
library(ggpubr)

setwd("~/Dropbox/GitHub/groundfish-o2-thresholds")

#Load functions
source("code/helper_funs.R")

#Load data
dat <- list.files(path = "data/processed_data/fish2", pattern = ".rds", full.names=T) %>%
  map(readRDS) %>% 
  bind_rows()

##Clean up data for plotting

#Remove any rows with necessary data missing
dat <- dat %>%
  drop_na(depth, mi1, temperature_C, salinity_psu, X, Y, year)

#Remove weird depths
dat <- filter(dat, depth>0)

#Remove oxygen outliers
dat <- filter(dat, O2_umolkg<1500)

#Add survey type
dat$survey_type <- ifelse(dat$survey=="iphc", "IPHC longline", "bottom trawl survey")

#Remove aleutian islands
dat <- filter(dat, region!="ai")

#Remove NBS
dat <- filter(dat, region!="nbs")

#Remove duplicates
dat <- unique(dat)

#Region labels
#Set region
dat$region <- factor(dat$region, levels=c("ebs", "goa", "bc", "cc"))
labs <- c("Eastern Bering Sea", "Gulf of Alaska", "British Columbia", "California Current")
names(labs) <- c("ebs", "goa", "bc", "cc")

#Plot themes
theme_set(theme_bw(base_size = 17))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

##Combine IPHC and bottom trawl surveys, positive catch only, for each species
##Plot barplot with IPHC and bottom trawl
#all species
ggplot(filter(dat, (survey!="iphc" & catch_weight>0)|(survey=="iphc"&(cpue_weight>0|cpue_count>0))), aes(x=year, fill=region, pattern=survey_type))+
  #stat_count(aes(fill=region, pattern=survey_type))+
  facet_wrap("common_name", ncol=4, scales="free_y", labeller=labeller(common_name=label_wrap_gen(20)))+
  scale_x_continuous(breaks=c(2009,2016,2023), limits=c(2008,2023))+
  geom_bar_pattern(
    stat = "count",
    colour = "black",            # Outline color
    pattern_fill = "black",      # Pattern stripe color
    pattern_angle = 45,
    pattern_density = 0.2,
    pattern_spacing = 0.05,
    pattern_size = 0.1)+
  scale_pattern_manual(
    name="data type",
    values = c("none", "stripe"))+
  xlab("Year")+
  ylab("Number of Observations")+
  theme(legend.position="top",  panel.spacing=unit(0, "pt"),legend.justification="center", legend.box.spacing = unit(0, "pt"))+
  guides(fill = guide_legend(nrow = 2), pattern=guide_legend(nrow=2,override.aes = list(pattern = c("none", "stripe"))))+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs)

ggsave(
  paste("output/plots/dat_availability/barplot_positiveonly_all.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 10,
  height = 13,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##All observations, example
example <- filter(dat, common_name=="sablefish")
#Labels
example$survey <- factor(example$survey, levels=c("iphc", "afsc_bsai", "afsc_goa" ,"dfo", "nwfsc")) 
labs2 <- c("IPHC Longline", "NOAA BS Bottom Trawl", "NOAA GOA Bottom Trawl", "DFO BC Bottom Trawl", "NOAA WC Bottom Trawl")
names(labs2) <- c("iphc", "afsc_bsai", "afsc_goa" ,"dfo", "nwfsc")

theme_set(theme_bw(base_size = 20))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

plot1 <- ggplot(example, aes(x=year, fill=survey))+
  stat_count(aes(fill=survey))+
  facet_wrap("region", ncol=2, scales="free_y", labeller=labeller(region=labs))+
  scale_x_continuous(breaks=c(2009,2016,2023), limits=c(2008,2023))+
  scale_pattern_manual(
    name="data type",
    values = c("none", "stripe"))+
  xlab("Year")+
  ylab("Number of Observations")+
  theme(legend.position="top")+
  guides(fill = guide_legend(nrow = 2))+
  scale_fill_manual(values=c("#f2cc84","#B9D7D9", "#a7ba42", "darkslategrey", "#FAD2E1"), drop=FALSE, labels=labs2)

plot1

 ggsave(
  paste("output/plots/data_available_overall.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 8.5,
  height = 7,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

###Map of all data available
# setup up mapping ####
map_data <- rnaturalearth::ne_countries(scale = "large",
                                        returnclass = "sf",
                                        continent = "North America")

us_coast_proj <- sf::st_transform(map_data, crs = 32610)

#Plot map
plot2 <- ggplot(us_coast_proj) + geom_sf() +
 geom_point(filter(example, survey_type!="bottom trawl survey"), mapping=aes(x=X*1000, y=Y*1000,colour=survey), size=0.3, alpha=0.5)+
 geom_point(filter(example, survey_type=="bottom trawl survey"), mapping=aes(x=X*1000, y=Y*1000,colour=survey), size=0.3, alpha=1)+
  scale_x_continuous(breaks=c(-150,-135,-120), limits=c(min(example$X)*1000, max(example$X)*1000))+
  ylim(min(example$Y)*1000, max(example$Y)*1000)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position=c(0.4,0.3))+
  guides(color = guide_legend(override.aes = list(size=3, alpha=1)))+
  scale_color_manual(values=c( "#f2cc84","#B9D7D9", "#a7ba42", "darkslategrey", "#FAD2E1"), drop=FALSE, labels=labs2)

plot2

ggsave(
  paste("output/plots/data_available_map.png"),
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

#Combine
#Combine plots
ggarrange(plot1, plot2, common.legend=T)

ggsave(
  paste("output/plots/data_available_combined.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 15,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#Other option, with above and below and colored
summary <- example %>%
  group_by(region, survey_type, year) %>%
  summarize(count=n())
summary$count <- ifelse(summary$survey_type=="IPHC longline", summary$count*-1, summary$count)

theme_set(theme_bw(base_size = 30))
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank())

plot3 <- ggplot(summary, aes(x=year, y=count, fill=region))+
  geom_col()+
 # facet_wrap("region", ncol=2, scales="free_y", labeller=labeller(region=labs))+
  #scale_x_continuous(breaks=c(2009,2016,2023), limits=c(2008,2023))+
  scale_y_continuous(breaks=c(-1000,-500,0,500,1000),labels=c(1000,500,0,500,1000), limits=c(-1200,1200))+
  xlab("Year")+
  ylab("Number of Observations")+
  guides(fill = guide_legend(nrow = 2))+
  scale_fill_manual(values=c("#88CCEE", "#999933", "#44AA99","#CC6677"), drop=FALSE, labels=labs)+
  theme(legend.position="top", legend.title=element_blank())+
  geom_hline(yintercept=0, linetype="solid")+
  annotate("text", x = -Inf, y = Inf, label = paste("Bottom Trawl Data"), vjust = 2, hjust = -0.15, size=10)+
  annotate("text", x = -Inf, y = -Inf, label = paste("IPHC Longline Data"), vjust=-0.5, hjust=-0.15, size=10)

plot3
ggsave(
  paste("output/plots/data_available_splitbar.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 15,
  height = 8.5,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

#Combine second option
plot4 <- ggplot(us_coast_proj) + geom_sf() +
  geom_point(filter(example, survey_type!="bottom trawl survey"), mapping=aes(x=X*1000, y=Y*1000,colour=survey), size=0.3, alpha=0.5)+
  geom_point(filter(example, survey_type=="bottom trawl survey"), mapping=aes(x=X*1000, y=Y*1000,colour=survey), size=0.3, alpha=1)+
  scale_x_continuous(breaks=c(-150,-135,-120), limits=c(min(example$X)*1000, max(example$X)*1000))+
  ylim(min(example$Y)*1000, max(example$Y)*1000)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(legend.position=c("top"), legend.title=element_blank())+
  guides(color = guide_legend(override.aes = list(size=6, alpha=1), nrow=3))+
  scale_color_manual(values=c("#f2cc84","lightsteelblue", "seagreen", "darkslategrey", "#FAD2E1"), drop=FALSE, labels=labs2)

plot3+plot4+plot_annotation(tag_levels='A')

ggsave(
  paste("output/plots/data_available_combined2.png"),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = 22,
  height = 16,
  units = c("in"),
  dpi = 600,
  limitsize = TRUE, bg="white"
)

##Make a table
#Number of positive catches per species per region
table <- dat %>% filter(survey_type!="iphc" & catch_weight>0)%>% 
  group_by(common_name, region) %>% 
  summarize(N_pos = length(catch_weight))
table <- pivot_wider(table, names_from=region, id_cols=common_name, values_from=c(N_pos))
#Species table
species <- read_excel("data/species_table.xlsx")
species$common_name <- tolower(species$common_name)
#Combine
species <- left_join(species, table)

##IPHC data
table2 <- dat %>% filter(survey=="iphc"& (cpue_count>0|cpue_weight>0))%>% 
  mutate(cpue=ifelse(is.na(cpue_count), cpue_weight, cpue_count))%>% 
  group_by(common_name, region) %>% 
  summarize(N_pos_IPHC = length(cpue))
  
table2 <- pivot_wider(table2, names_from=region, id_cols=common_name, values_from=c(N_pos_IPHC))
colnames(table2) <- c("common_name", "iphc_cc", "iphc_bc", "iphc_goa", "iphc_ebs")

species <- left_join(species, table2)

#save
write.csv(species, file="output/species_summary.csv")

##Models fit to data
#Plot
list <- unique(dat$common_name)
dat_names <- c("cc", "bc", "goa", "ebs", "coastwide", "cc _iphc", "bc _iphc", "goa _iphc", "ebs _iphc", "coastwide _iphc")
output_folder <- "region_comp"

for(i in 1:length(list)) {
  this_species = list[i]
  print(this_species)
  for(h in 1:length(dat_names)){
    this_dat <- dat_names[h]
    print(this_dat)
    dat2plot <- try(readRDS(file = paste0("output/",output_folder, "/", this_species, "_", this_dat, "_dat.rds")))
    if(is.data.frame(dat2plot)){
      #Labels
      dat2plot$survey <- factor(dat2plot$survey, levels=c("iphc", "afsc_bsai", "afsc_goa" ,"dfo", "nwfsc")) 
      labs2 <- c("IPHC Longline", "NOAA Bering Sea Bottom Trawl", "NOAA Gulf of Alaska Bottom Trawl", "DFO British Columbia Bottom Trawl", "NOAA West Coast Bottom Trawl")
      names(labs2) <- c("iphc", "afsc_bsai", "afsc_goa" ,"dfo", "nwfsc")
      ggplot(us_coast_proj) + geom_sf() +
        geom_point(filter(dat2plot,catch>0),mapping=aes(x=X*1000, y=Y*1000,colour=survey), size=0.1)+
        scale_x_continuous(breaks=c(-145,-120), limits=c(min(dat2plot$X)*1000, max(dat2plot$X)*1000))+        
        ylim(min(dat2plot$Y)*1000, max(dat2plot$Y)*1000)+
        facet_wrap("year", ncol=5)+
        theme_minimal(base_size=15)+
        xlab("Longitude")+
        ylab("Latitude")+
        theme(legend.position="top")+
        ggtitle(paste(this_species, this_dat, sep=" "))+
        scale_color_manual(values=c("#f2cc84","lightsteelblue", "seagreen", "darkslategrey", "#FAD2E1"), drop=FALSE, labels=labs2)+
        guides(color=guide_legend(override.aes=list(alpha=1, size=4), nrow=3))
      
      
      ggsave(
        paste0(output_folder, "/", "plots/data_fit_mapping/map_",this_species,"_", this_dat,".png"),
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
  }
}

