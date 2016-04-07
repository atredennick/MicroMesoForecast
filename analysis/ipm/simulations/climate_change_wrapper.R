#### Wrapper for running IPM climate change scenarios by species
####
#### Andrew Tredennick
#### 2-11-2014
#### Updated: 4-6-2016


do_change <- 1 # Determines the level of climate change programmatically


####
####  PRELIMINARIES
####
##  Use mean parameter values?
use_mean_params <- TRUE

##  Set working directory programmatically
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos") # modify as needed
setwd(paste(root,"/MicroMesoForecast/analysis/ipm/simulations",sep="")) # modify as needed 

##  Directory paths
if(use_mean_params==TRUE){
  resultspath <- "./results/climchange_meanparams/"
}
if(use_mean_params==FALSE){
  resultspath <- "./results/climchange_varyparams/"
}

##  Global parameters for the simulations
tlimit <- 250  # number of years to simulate
burn.in <- 50    # years to cut before calculations
spp_list <- c("BOGR","HECO","PASM","POSE")
n_spp <- length(spp_list)

##  Degree of climate change (can be a vector)
# clim_change <- c(0.01, 0.1, 0.2, 0.3) # climate change percentage
clim_change <- 0.01 # climate change percentage in decimal form
perc_change <- clim_change*100 # climate change percentage in %age form
clim_alt <- clim_change[do_change] # can be set programmatically if necessary
filetag <- paste(perc_change[do_change], ".csv", sep="")

##  Load climate scalers
growth_clim_scalers <- readRDS("../../growth_all_clim_scalers.RDS")
surv_clim_scalers <- readRDS("../../survival_all_clim_scalers.RDS")
rec_clim_scalers <- readRDS("../../recruitment_all_clim_scalers.RDS")

## Load climate and random year effect sequences
yrSave <- readRDS("../../random_year_effects_sequence.rds")
climYr <- readRDS("../../climate_year_sequence.rds")

## Fetch calendar years
raw_data <- read.csv("../../speciesData/quadAllCover.csv")
years <- unique(raw_data$year)
years <- years[1:(length(years)-1)] #lop off 1945 since no climate for that year

##  Fetch alpha values (needed to calculate neighborhood crowding)
##  Alpha values come directly from Chu and Adler (2015, Ecological Monographs)
alpha_grow <- read.csv("../../alpha_list_growth.csv")
alpha_surv <- read.csv("../../alpha_list_survival.csv")



####
####  SOURCE FUNCTION FILES FOR VITAL RATES AND PARAMETER IMPORT
####
source("vital_rate_ipm_functions.R", echo = FALSE) # for IPM

##  Functions for fetching vital rate parameters
if(use_mean_params==TRUE){
  source("../vitalRateRegs/survival/import2ipm_means.R", echo = FALSE)
  source("../vitalRateRegs/growth/import2ipm_means.R", echo = FALSE)
  source("../vitalRateRegs/recruitment/import2ipm_means.R", echo = FALSE)
}
if(use_mean_params==FALSE){
  source("../vitalRateRegs/survival/import2ipm.R", echo = FALSE)
  source("../vitalRateRegs/growth/import2ipm.R", echo = FALSE)
  source("../vitalRateRegs/recruitment/import2ipm.R", echo = FALSE) 
}



####
####  LOAD CLIMATE DATA (UNSTANDARDIZED)
####
clim_data_obs <- read.csv("../../weather/Climate.csv")
# Subset and reorder to match regression param import
clim_data_obs <- clim_data_obs[,c("year", "pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")]



####
####  FUNCTION FOR CLIMATE CHANGE IPM SIMULATIONS
####
climsimIPM <- function(clim_alt, clim_var, clim_data_obs, do_spp){
  doSpp <- spp_list[do_spp]
  
  if(is.na(clim_var)==TRUE){
    clim_data <- clim_data_obs
  }
  
  if(is.na(clim_var)==FALSE){
    # Perturb climate time series
    clim_data <- clim_data_obs # re-assigns the temporary clim_data object
    for(ivar in 1:length(clim_var)){
      vars <- grep(clim_var[ivar],names(clim_data_obs)) # indices for ppt columns
      tmp1 <- clim_alt*colMeans(clim_data_obs) # vector of ppt additions
      # matrix of additions to match DIMS of observed climate
      tmp1 <- matrix(tmp1,NROW(clim_data_obs),NCOL(clim_data_obs),byrow=T) 
      if(colMeans(tmp1)[1]!=tmp1[1,1]) { stop("climate change matrix mis-aligned") }
      clim_data[,vars] <- clim_data_obs[,vars]+tmp1[,vars] # applies ppt additions
    }
    # Make sure climate perturbations have been applied
    if(mean(clim_data$pptLag)==mean(clim_data_obs$pptLag)) {
      stop("climate perturbation not applied")
    }
  }
  
  # Retriev spp-specific scalers
  growth_clim_scalers <- subset(growth_clim_scalers, yearout==1 & species==doSpp)
  surv_clim_scalers <- subset(surv_clim_scalers, yearout==1 & species==doSpp)
  rec_clim_scalers <- subset(rec_clim_scalers, yearout==1 & species==doSpp)
  Gscalers <- growth_clim_scalers[,c("means","sds")]
  Sscalers <- surv_clim_scalers[,c("means","sds")]
  Rscalers <- rec_clim_scalers[,c("means","sds")]
  
  # Run the IPM
  source("ipm_climate_simulations.R", echo = FALSE)
  
  # Throw error if something weird happens
  if(min(covSave)<0.00001) { stop("really low values; check IPM") }
  if(max(covSave)>1) { stop("cover greater than 1; check IPM") }
  
  return(covSave)
}


####
#### RUN CLIMATE CHANGE SIMULATIONS
####
## Simulation 1: Observed Climate
clim_var <- NA
for(ispp in 1:length(spp_list)){
  covSave <- climsimIPM(clim_alt, clim_var, clim_data_obs, do_spp = ispp)
  
  # Set output filename
  filenameid <- paste(clim_var, collapse="")
  outfile <- paste0(resultspath,doSpp,"_ipm_cover_", filenameid, "Change.csv")
  climtag <- clim_var # tag for current climate simulations
  
  # Save cover time series (covSave from IPM)
  output <- data.frame("time"=burn.in:tlimit,"cover"=covSave[burn.in:tlimit])
  output$species <- doSpp
  output$climsim <- filenameid
  # write.table(output,outfile,row.names=F,sep=",")
  cat(paste("DONE WITH", spp_list[ispp], "FOR OBSERVED CLIMATE SIMULATION"))
  cat("\n") # makes a space before next line printing
}


####
####  SIMULATION 2: INCREASE IN PRECIPITATION
####
# Perturb climate time series
vars <- grep("ppt",names(clim_data_obs)) # indices for ppt columns
tmp1 <- clim_alt*colMeans(clim_data_obs) # vector of ppt additions

# matrix of additions to match DIMS of observed climate
tmp1 <- matrix(tmp1,NROW(clim_data_obs),NCOL(clim_data_obs),byrow=F) 

clim_data <- clim_data_obs # re-assigns the temporary clim_data object
clim_data[,vars] <- clim_data_obs[,vars]+tmp1[,vars] # applies ppt additions

for(ss in 1:n_spp){
  doSpp <- spp_list[ss] # set the species
  
  # Retriev spp-specific scalers
  growth_clim_scalers <- subset(growth_clim_scalers, yearout==1 & species==doSpp)
  surv_clim_scalers <- subset(surv_clim_scalers, yearout==1 & species==doSpp)
  rec_clim_scalers <- subset(rec_clim_scalers, yearout==1 & species==doSpp)
  Gscalers <- growth_clim_scalers[,c("means","sds")]
  Sscalers <- surv_clim_scalers[,c("means","sds")]
  Rscalers <- rec_clim_scalers[,c("means","sds")]
  
  # Set output filenames
  outfile <- paste0(resultspath,doSpp,"_ipm_cover_pptChange.csv")
  climtag <- "ppt" # tag for current climate simulations
  
  # Make sure climate perturbations have been applied
  if(mean(clim_data$pptLag)==mean(clim_data_obs$pptLag)) {
    stop("climate perturbation not applied")
  }
  
  # Run the IPM
  source("ipm_climate_simulations.R", echo = FALSE)
  
  # Throw error if something weird happens
  if(min(covSave)<0.00001) { stop("really low values; check IPM") }
  if(max(covSave)>1) { stop("cover greater than 1; check IPM") }
  
  # Save cover time series (covSave from IPM)
  output <- data.frame("time"=burn.in:tlimit,"cover"=covSave[burn.in:tlimit])
  output$species <- doSpp
  output$climsim <- climtag
  write.table(output,outfile,row.names=F,sep=",")
  cat(paste("DONE WITH", doSpp, "FOR INCREASED PRECIPITATION SIMULATION"))
  cat("\n") # makes a space before next line printing
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
  tmp1 <- clim_alt*colMeans(clim_data)
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
    outfile1<-paste("./results/climchange_varyparams/",doSpp,"_ipm_cover_tempChange", filetag ,sep="")
    outfile2<-paste("./results/climchange_varyparams/",doSpp,"_ipm_density_tempChange", filetag ,sep="")
    outfile3<-paste("./results/climchange_varyparams/",doSpp,"_ipm_stableSize_tempChange", filetag ,sep="")
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
  tmp1 <- clim_alt*colMeans(clim_data)
  tmp1 <- matrix(tmp1,NROW(clim_data),NCOL(clim_data),byrow=T)
  clim_data[,vars]=clim_data[,vars]+tmp1[,vars]
  # and again
  vars <- grep("ppt",names(clim_data))
  tmp1 <- clim_alt*colMeans(clim_data)
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
    outfile1<-paste("./results/climchange_varyparams/",doSpp,"_ipm_cover_temppptChange", filetag ,sep="")
    outfile2<-paste("./results/climchange_varyparams/",doSpp,"_ipm_density_temppptChange", filetag ,sep="")
    outfile3<-paste("./results/climchange_varyparams/",doSpp,"_ipm_stableSize_temppptChange", filetag ,sep="")
    climtag <- "zppttemp"
    source("ipm_climate_simulations.R", echo = FALSE)
    print(paste("DONE WITH", doSpp, "FOR TEMPERATURE+PRECIPITATION CHANGE SIMULATION"))
  }
} # end climate change wrapper

