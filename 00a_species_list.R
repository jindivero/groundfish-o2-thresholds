library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggplot2)
library(tidync)
library(sdmTMB)
library(marmap)
library(sf)

setwd("~/Dropbox/GitHub/groundfish-o2-sdm")
source("code/helper_funs.R")
source("code/util_funs.R")

##Edit species list
#load, clean, and join data
catch <- readRDS("data/fish_raw/BC/pbs-catch.rds")

itis <- unique(catch$itis)
itis <- itis[!is.na(itis)]
scientific_names <-  taxize::classification(itis, db="itis")
itis <- as.data.frame(itis)
itis <- drop_na(itis)

sci_names <- list()
for (i in 1:nrow(itis)){
  x <- as.data.frame(as.matrix(scientific_names[[i]]))
  last_row <- as.data.frame(x[nrow(x), ])
  name <- last_row$name
  name <-tolower(name)
  sci_names[[i]] <- name
}

names <- as.data.frame(unlist(sci_names))
itis$scientific_name <- names[,1]

#Add common name
species <- catch[,c("itis", "species_common_name")]
species <- unique(species)
species <- drop_na(species)
itis <- left_join(itis, species)
itis$common_name <- itis$species_common_name
itis$species_common_name <- NULL

#Rename dogfish
itis$common_name <- ifelse(itis$common_name=="north pacific spiny dogfish", "spiny dogfish", itis$common_name)

#Save

saveRDS(itis, file="data/fish_raw/BC/species-table.rds")
