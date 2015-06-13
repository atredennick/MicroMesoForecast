###############################################################
#### Wrapper for running climate change scenarios by species
####
#### Andrew Tredennick
#### 2-11-2014
####

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
setwd(paste(root,"/MicroMesoForecast/analysis/ipm/simulations",sep="")); # modify as needed 

# Set some global parameters for the simulations
tlimit<-2500  ## number of years to simulate
burn.in<-500    # years to cut before calculations
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)

# Get climate and random year effect sequences
yrSave <- readRDS("../../random_year_effects_sequence.rds")
climYr <- readRDS("../../climate_year_sequence.rds")
# paramSeq <- readRDS("../../random_chain_iter_sequence.RDS")

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
#### Simulation 1. Observed climate sequence -------------------------------------
####
# Set climate for observed climate run
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
clim_data[2:6] <- scale(clim_data[2:6], center = TRUE, scale = TRUE) # standardize

# Loop through species
# n_spp <- 1 #for test
for(ss in 1:n_spp){
  doSpp <- spp_list[ss]
  outfile1<-paste("./results/climchange_varyparams/",doSpp,"_ipm_cover_obsClimate.csv",sep="")
  outfile2<-paste("./results/climchange_varyparams/",doSpp,"_ipm_density_obsClimate.csv",sep="")
  outfile3<-paste("./results/climchange_varyparams/",doSpp,"_ipm_stableSize_obsClimate.csv",sep="")
  climtag <- "obs"
  source("ipm_climate_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR OBSERVED CLIMATE SIMULATION"))
}



####
#### Simulation 2. 1% increase in precipitation -------------------------------------
####
# Set climate for 1% ppt increase run
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
# clim_data[2:5] <- scale(clim_data[2:5], center = TRUE, scale = TRUE) # standardize
clim_avg <- apply(X = clim_data, MARGIN = 2, FUN = mean)
clim_sd <- apply(X = clim_data, MARGIN = 2, FUN = sd)

# Get just the variables of interest
vars <- grep("ppt",names(clim_data))
tmp1 <- 0.01*colMeans(clim_data)
tmp1 <- matrix(tmp1,NROW(clim_data),NCOL(clim_data),byrow=T)
clim_data[,vars]=clim_data[,vars]+tmp1[,vars]

# Now scale based on perturbed or regular data, depending on scenario
clim_data["pptLag"] <- (clim_data["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
clim_data["ppt1"] <- (clim_data["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
clim_data["ppt2"] <- (clim_data["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
clim_data["TmeanSpr1"] <- (clim_data["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
clim_data["TmeanSpr2"] <- (clim_data["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]

# Loop through species
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)
# n_spp <- 1
for(ss in 1:n_spp){
  doSpp <- spp_list[ss]
  outfile1<-paste("./results/climchange_varyparams/",doSpp,"_ipm_cover_pptChange.csv",sep="")
  outfile2<-paste("./results/climchange_varyparams/",doSpp,"_ipm_density_pptChange.csv",sep="")
  outfile3<-paste("./results/climchange_varyparams/",doSpp,"_ipm_stableSize_pptChange.csv",sep="")
  climtag <- "ppt"
  source("ipm_climate_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR PRECIPITATION CHANGE SIMULATION"))
}



####
#### Simulation 3. 1% increase in temperature -------------------------------------
####
# Set climate for 1% temp increase run
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year", "pptLag", "ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
# clim_data[2:5] <- scale(clim_data[2:5], center = TRUE, scale = TRUE) # standardize
clim_avg <- apply(X = clim_data, MARGIN = 2, FUN = mean)
clim_sd <- apply(X = clim_data, MARGIN = 2, FUN = sd)

# Get just the variables of interest
vars <- grep("Tmean",names(clim_data))
tmp1 <- 0.01*colMeans(clim_data)
tmp1 <- matrix(tmp1,NROW(clim_data),NCOL(clim_data),byrow=T)
clim_data[,vars]=clim_data[,vars]+tmp1[,vars]

# Now scale based on perturbed or regular data, depending on scenario
clim_data["pptLag"] <- (clim_data["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
clim_data["ppt1"] <- (clim_data["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
clim_data["ppt2"] <- (clim_data["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
clim_data["TmeanSpr1"] <- (clim_data["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
clim_data["TmeanSpr2"] <- (clim_data["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]

# Loop through species
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)
# n_spp <- 1
for(ss in 1:n_spp){
  doSpp <- spp_list[ss]
  outfile1<-paste("./results/climchange_varyparams/",doSpp,"_ipm_cover_tempChange.csv",sep="")
  outfile2<-paste("./results/climchange_varyparams/",doSpp,"_ipm_density_tempChange.csv",sep="")
  outfile3<-paste("./results/climchange_varyparams/",doSpp,"_ipm_stableSize_tempChange.csv",sep="")
  climtag <- "temp"
  source("ipm_climate_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR TEMPERATURE CHANGE SIMULATION"))
}



####
#### Simulation 4. 1% increase in temperature and precipitation -------------------------------------
####
# Set climate for 1% temp & ppt increase run
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year", "pptLag", "ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
# clim_data[2:5] <- scale(clim_data[2:5], center = TRUE, scale = TRUE) # standardize
clim_avg <- apply(X = clim_data, MARGIN = 2, FUN = mean)
clim_sd <- apply(X = clim_data, MARGIN = 2, FUN = sd)

# Get just the variables of interest
vars <- grep("Tmean",names(clim_data))
tmp1 <- 0.01*colMeans(clim_data)
tmp1 <- matrix(tmp1,NROW(clim_data),NCOL(clim_data),byrow=T)
clim_data[,vars]=clim_data[,vars]+tmp1[,vars]
# and again
vars <- grep("ppt",names(clim_data))
tmp1 <- 0.01*colMeans(clim_data)
tmp1 <- matrix(tmp1,NROW(clim_data),NCOL(clim_data),byrow=T)
clim_data[,vars]=clim_data[,vars]+tmp1[,vars]

# Now scale based on perturbed or regular data, depending on scenario
clim_data["pptLag"] <- (clim_data["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
clim_data["ppt1"] <- (clim_data["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
clim_data["ppt2"] <- (clim_data["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
clim_data["TmeanSpr1"] <- (clim_data["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
clim_data["TmeanSpr2"] <- (clim_data["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]

# Loop through species
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)
# n_spp <- 1
for(ss in 1:n_spp){
  doSpp <- spp_list[ss]
  outfile1<-paste("./results/climchange_varyparams/",doSpp,"_ipm_cover_temppptChange.csv",sep="")
  outfile2<-paste("./results/climchange_varyparams/",doSpp,"_ipm_density_temppptChange.csv",sep="")
  outfile3<-paste("./results/climchange_varyparams/",doSpp,"_ipm_stableSize_temppptChange.csv",sep="")
  climtag <- "zppttemp"
  source("ipm_climate_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR TEMPERATURE+PRECIPITATION CHANGE SIMULATION"))
}

