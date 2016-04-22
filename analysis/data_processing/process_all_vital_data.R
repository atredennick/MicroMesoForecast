##  R Script to Process and Transform Raw Data into Dataframes Used
##  for Statistical Modeling and Model Simulations.
##
##  Author: Andrew Tredennick
##  Last update: 4-21-2016

# Clear the workspace
rm(list=ls())


####
####  LIBRARIES
####
library("plyr")
library("reshape2")



####
####  PRELIMINARIES
####
species <- c("BOGR", "HECO", "PASM", "POSE")
path2data <- "./speciesData/"
outpath <- "../processed_data/"
surv_outfile <- "surv_with_weather.RDS"
grow_outfile <- "grow_with_weather.RDS"
rec_outfile <- "rec_with_weather.RDS"
cover_outfile <- "cover_with_weather.RDS"
weather <- read.csv("../weather/Climate.csv")
my_covariates <- c("year","pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
weather <- weather[,my_covariates]



####
####  SOURCE FILES
####
source("dataExtractQuad.R") # Adds in zeros for all measured quads; rms bad bogr quadyrs
source("format_quad_data_fxn.R") # Function for formatting quadrat data



####
####  QUADRAT PROPORTIONAL COVER DATA
####
quad_dat <- read.csv(paste0(path2data,"quadAllCover.csv"))
out_quad_dat <- format_quad_data(quad_dat, weather)
saveRDS(out_quad_dat, paste0(outpath,cover_outfile))



####
####  SURVIVAL DATA
####
surv_data <- list() # empty object to collate survival data frames
for(do_species in species){
  if(do_species=="BOGR") { getfile <- paste0(path2data,do_species,"/edited/survD.csv") }
  if(do_species!="BOGR") { getfile <- paste0(path2data,do_species,"/survD.csv") }
  tmpsurv <- read.csv(getfile)
  tmpsurv$species <- do_species
  surv_data <- rbind(surv_data, tmpsurv)
}

# Get rid of the 1900 for year in a weather df copy; to match survival df year
weather_copy <- weather
weather_copy$year <- weather_copy$year-1900

# Merge in weather data
surv_data <- merge(surv_data,weather_copy)

# Extract group information
surv_data$Group <- as.factor(substr(surv_data$quad,1,1))

# Read in previously estimated crowding indices
c1 <- read.csv("BOGRsurvCrowding.csv")[,2:3]
c1$species <- sppList[1]
c2 <- read.csv("HECOsurvCrowding.csv")[,2:3]
c2$species <- sppList[2]
c3 <- read.csv("PASMsurvCrowding.csv")[,2:3]
c3$species <- sppList[3]
c4 <- read.csv("POSEsurvCrowding.csv")[,2:3]
c4$species <- sppList[4]
crowd <- rbind(c1,c2,c3,c4)
colnames(crowd) <- c("W", "X", "species")

# Merge crowding and growth data
surv_data <- merge(surv_data, crowd, by=c("species", "X"))

# Save the survival data frame
saveRDS(surv_data, paste0(outpath,surv_outfile))






