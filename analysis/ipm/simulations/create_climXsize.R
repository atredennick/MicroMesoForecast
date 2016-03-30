##  R Script to Calculate All Possible SizeXClimate Effects for IPM Simulations
##  Uses observed weather data and species' size bins from IPM

# Clear The Workspace
rm(list=ls())


####
####  Load Libraries
####
library("reshape2")
library("plyr")



####
####  Set Directory Paths and Filenames
####
weather_file <- "../../weather/Climate.csv"
path2species <- "../../speciesData/"
path2scalers <- "../../"
scaler_files <- list.files(path2scalers)[grep("scalers", list.files(path2scalers))]
torms <- grep("qbm", scaler_files)
scaler_files <- scaler_files[-torms]
scaler_files <- scaler_files[c(1,3)]
vitals <- c("growth", "survival") # In order of scaler_files



####
####  Global Settings from IPM
####
spp_list<-c("BOGR","HECO","PASM","POSE")
Atotal=10000 #Area of 100cm x 100cm quadrat
bigM=c(75,50,10,50) #Set matrix dimension for each species
maxSize=c(2500,120,25,100) #Max size for each species

get_midpoints <- function(bigM_spp, maxSize_spp){
  # stuff for numerical approximation....
  # minimum (0.9*minimum size from data) and 
  # maximum sizes (1.1*maximum size from data)
  # log-transformation...
  # so that L, U, b, v, h are all log-transformed....
  L=log(0.2)
  U=log(maxSize_spp)*1.1     
  # b---boundary points. Note: b chops up the size interval (L-U) into bigM-equal-sized portions.
  b = L+c(0:bigM_spp)*(U-L)/bigM_spp #log-transformed
  # v---middle points (in Ellner's script, v is y); calculates the middle of each n-equal-sized portion.
  v = 0.5*(b[1:bigM_spp]+b[2:(bigM_spp+1)]) #log-transformed
  return(v)
}



####
####  Loop Over Species and Vital Rates to Make Look Up Tables
####
clim_data <- read.csv(weather_file)
years <- clim_data$year
all_spp_list <- list()
for(do_spp in spp_list){
  pathspp <- paste0(path2species,do_spp,"/")
  ii <- which(spp_list==do_spp)
  u <- get_midpoints(bigM_spp = bigM[ii], maxSize_spp = maxSize[ii])
  
  vital_list <- list()
  for(do_vital in vitals){
    vital_scalers <- readRDS(paste0(path2scalers,scaler_files[grep(do_vital, scaler_files)]))
    vital_scalers <- subset(vital_scalers, species==do_spp & yearout==1)
    # Check covariate scaling lengths
    if(do_vital!="recruitment" & nrow(vital_scalers)!=12) { 
      stop("check vital scalers data frame") }
    
    for(do_year in years){
      weather <- subset(clim_data, year==do_year)[,2:6]
      weather <- weather[,c("pptLag", "ppt1","ppt2","TmeanSpr1","TmeanSpr2")]
      weather$inter1 <- weather$ppt1*weather$TmeanSpr1
      weather$inter2 <- weather$ppt2*weather$TmeanSpr2
      
      # clim_mains <- (weather - vital_scalers[1:length(weather),"means"])/vital_scalers[1:length(weather),"sds"]
      clim_size <- as.data.frame(as.matrix(apply(expand.grid(as.numeric(weather[1:5]),u), 1, prod)))
      clim_size$means <- rep(vital_scalers[8:nrow(vital_scalers),"means"], times=length(u))
      clim_size$sds <- rep(vital_scalers[8:nrow(vital_scalers),"sds"], times=length(u))
      clim_size$climXsize <- with(clim_size, (V1-means)/sds)
      clim_size$midpoint <- rep(u, each=length(8:nrow(vital_scalers)))
      clim_size$covariate <- rep(c("pptLagSize", "ppt1Size", "ppt2Size", "TmeanSpr1Size", "TmeanSpr2Size"), times=length(u))
      colnames(clim_size)[1] <- "raw_cov_value"
      vital_list[[do_vital]][[as.character(do_year)]] <- clim_size
    } # End year loop
  } # End vital rate loop
  all_spp_list[[do_spp]] <- vital_list
} # End species loop

saveRDS(all_spp_list, "spp_year_climateXsize_scaled.RDS")

