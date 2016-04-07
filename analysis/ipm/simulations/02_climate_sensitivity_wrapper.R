###############################################################
#### Wrapper for running climate change scenarios by species
####
#### Andrew Tredennick
#### 2-11-2014
####

root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
setwd(paste(root,"/MicroMesoForecast/analysis/ipm/simulations",sep="")); # modify as needed 

####
##  Set up climate datasets
####
#Increase precipitation by 1%
clim_ppt <- read.csv("../../weather/Climate.csv")
clim_ppt <- clim_ppt[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
# clim_data[2:5] <- scale(clim_data[2:5], center = TRUE, scale = TRUE) # standardize
clim_avg <- apply(X = clim_ppt, MARGIN = 2, FUN = mean)
clim_sd <- apply(X = clim_ppt, MARGIN = 2, FUN = sd)
# Get just the variables of interest
vars <- grep("ppt",names(clim_ppt))
tmp1 <- 0.01*colMeans(clim_ppt)
tmp1 <- matrix(tmp1,NROW(clim_ppt),NCOL(clim_ppt),byrow=T)
clim_ppt[,vars]=clim_ppt[,vars]+tmp1[,vars]
# Now scale based on perturbed or regular data, depending on scenario
clim_ppt["pptLag"] <- (clim_ppt["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
clim_ppt["ppt1"] <- (clim_ppt["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
clim_ppt["ppt2"] <- (clim_ppt["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
clim_ppt["TmeanSpr1"] <- (clim_ppt["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
clim_ppt["TmeanSpr2"] <- (clim_ppt["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]


# Set climate for 1% temp increase run
clim_temp <- read.csv("../../weather/Climate.csv")
clim_temp <- clim_temp[,c("year", "pptLag", "ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
# clim_temp[2:5] <- scale(clim_temp[2:5], center = TRUE, scale = TRUE) # standardize
clim_avg <- apply(X = clim_temp, MARGIN = 2, FUN = mean)
clim_sd <- apply(X = clim_temp, MARGIN = 2, FUN = sd)
# Get just the variables of interest
vars <- grep("Tmean",names(clim_temp))
tmp1 <- 0.01*colMeans(clim_temp)
tmp1 <- matrix(tmp1,NROW(clim_temp),NCOL(clim_temp),byrow=T)
clim_temp[,vars]=clim_temp[,vars]+tmp1[,vars]
# Now scale based on perturbed or regular data, depending on scenario
clim_temp["pptLag"] <- (clim_temp["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
clim_temp["ppt1"] <- (clim_temp["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
clim_temp["ppt2"] <- (clim_temp["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
clim_temp["TmeanSpr1"] <- (clim_temp["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
clim_temp["TmeanSpr2"] <- (clim_temp["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]

clim_both <- read.csv("../../weather/Climate.csv")
clim_both <- clim_both[,c("year", "pptLag", "ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
# clim_both[2:5] <- scale(clim_both[2:5], center = TRUE, scale = TRUE) # standardize
clim_avg <- apply(X = clim_both, MARGIN = 2, FUN = mean)
clim_sd <- apply(X = clim_both, MARGIN = 2, FUN = sd)
# Get just the variables of interest
vars <- grep("Tmean",names(clim_both))
tmp1 <- 0.01*colMeans(clim_both)
tmp1 <- matrix(tmp1,NROW(clim_both),NCOL(clim_both),byrow=T)
clim_both[,vars]=clim_both[,vars]+tmp1[,vars]
# and again
vars <- grep("ppt",names(clim_both))
tmp1 <- 0.01*colMeans(clim_both)
tmp1 <- matrix(tmp1,NROW(clim_both),NCOL(clim_both),byrow=T)
clim_both[,vars]=clim_both[,vars]+tmp1[,vars]
# Now scale based on perturbed or regular data, depending on scenario
clim_both["pptLag"] <- (clim_both["pptLag"] - clim_avg["pptLag"])/clim_sd["pptLag"]
clim_both["ppt1"] <- (clim_both["ppt1"] - clim_avg["ppt1"])/clim_sd["ppt1"]
clim_both["ppt2"] <- (clim_both["ppt2"] - clim_avg["ppt2"])/clim_sd["ppt2"]
clim_both["TmeanSpr1"] <- (clim_both["TmeanSpr1"] - clim_avg["TmeanSpr1"])/clim_sd["TmeanSpr1"]
clim_both["TmeanSpr2"] <- (clim_both["TmeanSpr2"] - clim_avg["TmeanSpr2"])/clim_sd["TmeanSpr2"]


# Set some global parameters for the simulations
tlimit<-2500  ## number of years to simulate
burn.in<-500    # years to cut before calculations
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)

# Source this to generate common year sequences (they are read in below)
# Requires tlimit to be set
# source("get_year_sequence.R", echo = FALSE)

# Get climate and random year effect sequences
yrSave <- readRDS("../../random_year_effects_sequence.rds")
climYr <- readRDS("../../climate_year_sequence.rds")

# get calendar years
raw_data <- read.csv("../../speciesData/quadAllCover.csv")
years <- unique(raw_data$year)
years <- years[1:(length(years)-1)] #lop off 1945 since no climate for that year

# get alphas values (needed to calculate neighborhood crowding)
alpha_grow <- read.csv("../../alpha_list_growth.csv")
alpha_surv <- read.csv("../../alpha_list_survival.csv")

# Source functions
source("vital_rate_ipm_functions.R", echo = FALSE)
source("../vitalRateRegs/survival/import2ipm_means.R", echo = FALSE)
source("../vitalRateRegs/growth/import2ipm_means.R", echo = FALSE)
source("../vitalRateRegs/recruitment/import2ipm_means.R", echo = FALSE)
# source("../vitalRateRegs/survival/import2ipm.R", echo = FALSE)
# source("../vitalRateRegs/growth/import2ipm.R", echo = FALSE)
# source("../vitalRateRegs/recruitment/import2ipm.R", echo = FALSE)


####
##  Simulation 1. All climate unperturbed with all vital rates ----------------
####
# Set climate for observed climate run
clim_data <- read.csv("../../weather/Climate.csv")
clim_data <- clim_data[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
clim_data[2:6] <- scale(clim_data[2:6], center = TRUE, scale = TRUE) # standardize

# Loop through species
# n_spp <- 1 #for test
for(ss in 1:n_spp){
  doSpp <- spp_list[ss]
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_Baseline.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_Baseline.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_Baseline.csv",sep="")
  doPpt <- "none"
  doTemp <- "none"
  doBoth <- "none"
  simcode <- "all"
  vitalcode <- "all"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR CLIMATE SIMULATION"))
}



####
##  Simulation 2. Growth perturb ppt increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_GrowPpt.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_GrowPpt.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_GrowPpt.csv",sep="")
  doPpt <- "growth"
  doTemp <- "none"
  simcode <- "ppt"
  vitalcode <- "growth"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR GROW PRECIPITATION SIMULATION"))
}


####
##  Simulation 3. Growth perturb temp increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_GrowTemp.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_GrowTemp.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_GrowTemp.csv",sep="")
  doPpt <- "none"
  doTemp <- "growth"
  simcode <- "temp"
  vitalcode <- "growth"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR GROW TEMPERATURE SIMULATION"))
}


####
##  Simulation 4. Survival ppt increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_SurvPpt.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_SurvPpt.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_SurvPpt.csv",sep="")
  doPpt <- "survival"
  doTemp <- "none"
  simcode <- "ppt"
  vitalcode <- "survival"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR SURVIVAL PRECIPITATION SIMULATION"))
}


####
##  Simulation 5. Survival temp increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_SurvTemp.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_SurvTemp.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_SurvTemp.csv",sep="")
  doPpt <- "none"
  doTemp <- "survival"
  simcode <- "temp"
  vitalcode <- "survival"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR NO SURVIVAL TEMPERATURE SIMULATION"))
}


####
##  Simulation 6. Recruitment ppt increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_RecPpt.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_RecPpt.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_RecPpt.csv",sep="")
  doPpt <- "recruitment"
  doTemp <- "none"
  simcode <- "ppt"
  vitalcode <- "recruitment"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR RECRUITMENT PRECIPITATION SIMULATION"))
}


####
##  Simulation 7. Recruitment temp increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_RecTemp.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_RecTemp.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_RecTemp.csv",sep="")
  doPpt <- "none"
  doTemp <- "recruitment"
  simcode <- "temp"
  vitalcode <- "recruitment"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR RECRUITMENT TEMPERATURE SIMULATION"))
}


####
##  Simulation 8. Growth ppt+temp increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_GrowPptTemp.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_GrowPptTemp.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_GrowPptTemp.csv",sep="")
  doPpt <- "none"
  doTemp <- "none"
  doBoth <- "growth"
  simcode <- "ppttemp"
  vitalcode <- "growth"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR GROW PPT+TEMP SIMULATION"))
}

####
##  Simulation 9. Survival ppt+temp increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_SurvPptTemp.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_SurvPptTemp.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_SurvPptTemp.csv",sep="")
  doPpt <- "none"
  doTemp <- "none"
  doBoth <- "survival"
  simcode <- "ppttemp"
  vitalcode <- "survival"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR SURV PPT+TEMP SIMULATION"))
}

####
##  Simulation 10. Recruitment ppt+temp increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_RecPptTemp.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_RecPptTemp.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_RecPptTemp.csv",sep="")
  doPpt <- "none"
  doTemp <- "none"
  doBoth <- "recruitment"
  simcode <- "ppttemp"
  vitalcode <- "recruitment"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR REC PPT+TEMP SIMULATION"))
}


####
##  Simulation 11. GrowthSurv ppt increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_GrowSurvPpt.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_GrowSurvPpt.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_GrowSurvPpt.csv",sep="")
  doPpt <- "growth_surv"
  doTemp <- "none"
  simcode <- "ppt"
  vitalcode <- "growth_surv"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR REC PPT+TEMP SIMULATION"))
}

####
##  Simulation 12. GrowthRec ppt increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_GrowRecPpt.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_GrowRecPpt.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_GrowRecPpt.csv",sep="")
  doPpt <- "growth_rec"
  doTemp <- "none"
  simcode <- "ppt"
  vitalcode <- "growth_rec"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR REC PPT+TEMP SIMULATION"))
}

####
##  Simulation 13. SurvRec ppt increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_SurvRecPpt.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_SurvRecPpt.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_SurvRecPpt.csv",sep="")
  doPpt <- "surv_rec"
  doTemp <- "none"
  simcode <- "ppt"
  vitalcode <- "surv_rec"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR REC PPT+TEMP SIMULATION"))
}


####
##  Simulation 14. GrowthSurv temp increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_GrowSurvTemp.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_GrowSurvTemp.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_GrowSurvTemp.csv",sep="")
  doTemp <- "growth_surv"
  doPpt <- "none"
  simcode <- "temp"
  vitalcode <- "growth_surv"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR REC PPT+TEMP SIMULATION"))
}

####
##  Simulation 15. GrowthRec temp increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_GrowRecTemp.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_GrowRecTemp.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_GrowRecTemp.csv",sep="")
  doTemp <- "growth_rec"
  doPpt <- "none"
  simcode <- "temp"
  vitalcode <- "growth_rec"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR REC PPT+TEMP SIMULATION"))
}

####
##  Simulation 16. SurvRec temp increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_SurvRecTemp.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_SurvRecTemp.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_SurvRecTemp.csv",sep="")
  doTemp <- "surv_rec"
  doPpt <- "none"
  simcode <- "temp"
  vitalcode <- "surv_rec"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR REC PPT+TEMP SIMULATION"))
}



####
##  Simulation 17. GrowthSurv ppt+temp increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_GrowSurvPptTemp.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_GrowSurvPptTemp.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_GrowSurvPptTemp.csv",sep="")
  doPpt <- "none"
  doTemp <- "none"
  doBoth <- "growth_surv"
  simcode <- "ppttemp"
  vitalcode <- "growth_surv"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR REC PPT+TEMP SIMULATION"))
}

####
##  Simulation 18. GrowthRec ppt+temp increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_GrowRecPptTemp.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_GrowRecPptTemp.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_GrowRecPptTemp.csv",sep="")
  doPpt <- "none"
  doTemp <- "none"
  doBoth <- "growth_rec"
  simcode <- "ppttemp"
  vitalcode <- "growth_rec"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR REC PPT+TEMP SIMULATION"))
}

####
##  Simulation 19. SurvRec ppt+temp increase -------------------------------------
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
  outfile1<-paste("./results/climate_sensitivity/",doSpp,"_ipm_cover_SurvRecPptTemp.csv",sep="")
  outfile2<-paste("./results/climate_sensitivity/",doSpp,"_ipm_density_SurvRecPptTemp.csv",sep="")
  outfile3<-paste("./results/climate_sensitivity/",doSpp,"_ipm_stableSize_SurvRecPptTemp.csv",sep="")
  doPpt <- "none"
  doTemp <- "none"
  doBoth <- "surv_rec"
  simcode <- "ppttemp"
  vitalcode <- "surv_rec"
  source("ipm_climate_sensitivity_simulations.R", echo = FALSE)
  print(paste("DONE WITH", doSpp, "FOR REC PPT+TEMP SIMULATION"))
}




