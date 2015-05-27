###############################################################
#### Wrapper for running climate change scenarios by species
####
#### Andrew Tredennick
#### 2-11-2014
####

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
setwd(paste(root,"/MicroMesoForecast/montana/ipm/simulations",sep="")); # modify as needed 

# Set some global parameters for the simulations
tlimit<-1100  ## number of years to simulate
burn.in<-100    # years to cut before calculations
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)

# Source this to generate common year sequences (they are read in below)
# Requires tlimit to be set
source("get_year_sequence.R", echo = FALSE)

# Get climate and random year effect sequences
yrSave <- readRDS("random_year_effects_sequence.rds")
climYr <- readRDS("climate_year_sequence.rds")

# get calendar years
raw_data <- read.csv("../../speciesData/quadAllCover.csv")
years <- unique(raw_data$year)
years <- years[1:(length(years)-1)] #lop off 1945 since no climate for that year

# get alphas values (needed to calculate neighborhood crowding)
alpha_grow <- read.csv("../../alpha_list_growth.csv")
alpha_surv <- read.csv("../../alpha_list_survival.csv")

# Source functions
source("vital_rate_ipm_functions.R", echo = FALSE)
source("../vitalRateRegs/survival/import2ipm.R", echo = FALSE)
source("../vitalRateRegs/growth/import2ipm.R", echo = FALSE)
source("../vitalRateRegs/recruitment/import2ipm.R", echo = FALSE)


####
#### Simulation 1. All climate covariates -------------------------------------
####
# Set climate for observed climate run
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
clim_data[2:6] <- scale(clim_data[2:6], center = TRUE, scale = TRUE) # standardize

# Loop through species
# n_spp <- 1 #for test
for(ss in 1:n_spp){
  doSpp <- spp_list[ss]
  outfile1<-paste("./sensitivity_results/",doSpp,"_ipm_cover_allClimate.csv",sep="")
  outfile2<-paste("./sensitivity_results/",doSpp,"_ipm_density_allClimate.csv",sep="")
  outfile3<-paste("./sensitivity_results/",doSpp,"_ipm_stableSize_allClimate.csv",sep="")
  noPpt <- "none"
  noTemp <- "none"
  source("ipm_climate_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR ALL CLIMATE SIMULATION"))
}



####
#### Simulation 2. No growth precipitation coefficients -------------------------------------
####
# Set climate for observed climate run
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
clim_data[2:6] <- scale(clim_data[2:6], center = TRUE, scale = TRUE) # standardize

# Loop through species
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)
# n_spp <- 1
for(ss in 1:n_spp){
  doSpp <- spp_list[ss]
  outfile1<-paste("./sensitivity_results/",doSpp,"_ipm_cover_noGrowPpt.csv",sep="")
  outfile2<-paste("./sensitivity_results/",doSpp,"_ipm_density_noGrowPpt.csv",sep="")
  outfile3<-paste("./sensitivity_results/",doSpp,"_ipm_stableSize_noGrowPpt.csv",sep="")
  noPpt <- "growth"
  noTemp <- "none"
  source("ipm_climate_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR NO GROW PRECIPITATION SIMULATION"))
}


####
#### Simulation 3. No growth temperature coefficients -------------------------------------
####
# Set climate for observed climate run
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
clim_data[2:6] <- scale(clim_data[2:6], center = TRUE, scale = TRUE) # standardize

# Loop through species
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)
# n_spp <- 1
for(ss in 1:n_spp){
  doSpp <- spp_list[ss]
  outfile1<-paste("./sensitivity_results/",doSpp,"_ipm_cover_noGrowTemp.csv",sep="")
  outfile2<-paste("./sensitivity_results/",doSpp,"_ipm_density_noGrowTemp.csv",sep="")
  outfile3<-paste("./sensitivity_results/",doSpp,"_ipm_stableSize_noGrowTemp.csv",sep="")
  noPpt <- "none"
  noTemp <- "growth"
  source("ipm_climate_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR NO GROW TEMPERATURE SIMULATION"))
}


####
#### Simulation 4. No survival precipitation coefficients -------------------------------------
####
# Set climate for observed climate run
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
clim_data[2:6] <- scale(clim_data[2:6], center = TRUE, scale = TRUE) # standardize

# Loop through species
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)
# n_spp <- 1
for(ss in 1:n_spp){
  doSpp <- spp_list[ss]
  outfile1<-paste("./sensitivity_results/",doSpp,"_ipm_cover_noSurvPpt.csv",sep="")
  outfile2<-paste("./sensitivity_results/",doSpp,"_ipm_density_noSurvPpt.csv",sep="")
  outfile3<-paste("./sensitivity_results/",doSpp,"_ipm_stableSize_noSurvPpt.csv",sep="")
  noPpt <- "survival"
  noTemp <- "none"
  source("ipm_climate_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR NO SURVIVAL PRECIPITATION SIMULATION"))
}


####
#### Simulation 5. No survival temperature coefficients -------------------------------------
####
# Set climate for observed climate run
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
clim_data[2:6] <- scale(clim_data[2:6], center = TRUE, scale = TRUE) # standardize

# Loop through species
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)
# n_spp <- 1
for(ss in 1:n_spp){
  doSpp <- spp_list[ss]
  outfile1<-paste("./sensitivity_results/",doSpp,"_ipm_cover_noSurvTemp.csv",sep="")
  outfile2<-paste("./sensitivity_results/",doSpp,"_ipm_density_noSurvTemp.csv",sep="")
  outfile3<-paste("./sensitivity_results/",doSpp,"_ipm_stableSize_noSurvTemp.csv",sep="")
  noPpt <- "none"
  noTemp <- "survival"
  source("ipm_climate_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR NO SURVIVAL TEMPERATURE SIMULATION"))
}


####
#### Simulation 6. No recruitment precipitation coefficients -------------------------------------
####
# Set climate for observed climate run
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
clim_data[2:6] <- scale(clim_data[2:6], center = TRUE, scale = TRUE) # standardize

# Loop through species
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)
# n_spp <- 1
for(ss in 1:n_spp){
  doSpp <- spp_list[ss]
  outfile1<-paste("./sensitivity_results/",doSpp,"_ipm_cover_noRecPpt.csv",sep="")
  outfile2<-paste("./sensitivity_results/",doSpp,"_ipm_density_noRecPpt.csv",sep="")
  outfile3<-paste("./sensitivity_results/",doSpp,"_ipm_stableSize_noRecPpt.csv",sep="")
  noPpt <- "recruitment"
  noTemp <- "none"
  source("ipm_climate_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR NO RECRUITMENT PRECIPITATION SIMULATION"))
}


####
#### Simulation 7. No recruitment temperature coefficients -------------------------------------
####
# Set climate for observed climate run
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
clim_data[2:6] <- scale(clim_data[2:6], center = TRUE, scale = TRUE) # standardize

# Loop through species
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)
# n_spp <- 1
for(ss in 1:n_spp){
  doSpp <- spp_list[ss]
  outfile1<-paste("./sensitivity_results/",doSpp,"_ipm_cover_noRecTemp.csv",sep="")
  outfile2<-paste("./sensitivity_results/",doSpp,"_ipm_density_noRecTemp.csv",sep="")
  outfile3<-paste("./sensitivity_results/",doSpp,"_ipm_stableSize_noRecTemp.csv",sep="")
  noPpt <- "none"
  noTemp <- "recruitment"
  source("ipm_climate_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR NO RECRUITMENT TEMPERATURE SIMULATION"))
}




