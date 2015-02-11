#######################################################
#### Script to produce random sequence of years
#### for climate draws and random year effect draws


# single species IPM with climate covariates
# this script simulates equilibrium cover

#============================================================
# (I) INPUTS
#============================================================

#set working directory
root=ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
setwd(paste(root,"/MicroMesoForecast/montana/ipm/simulations",sep="")); # modify as needed 

# get calendar years
raw_data <- read.csv("../../speciesData/quadAllCover.csv")
years <- unique(raw_data$year)
years <- years[1:(length(years)-1)] #lop off 1945 since no climate for that year

# produce 2500 year time-series
tlimit <- 2500
yrSave <- sample(years,tlimit,replace=T) 
climYr <- sample(years,tlimit,replace=T) 

# save the time series
saveRDS(yrSave, "random_year_effects_sequence.rds")
saveRDS(climYr, "climate_year_sequence.rds")



