##  R Script for IPM forecasts of final observation year from all possible
##  lag years. IPM is initialized with cover at year t-x in quadrat q, and 
##  then projected using observed climate to year t. We can then assess the 
##  model's forecast horizon and compare it to the IBM.
##
##  Author: Andrew Tredennick
##  Last update: 4-19-2016


# Clear the workspace
rm(list=ls())



####
####  PRELIMINARIES
####
spp_list <- c("BOGR","HECO","PASM","POSE") # all Montana species
reps <- 50 # number of times to simulate each quadrat x year 

##  Load climate scalers
growth_clim_scalers <- readRDS("../../growth_all_clim_scalers.RDS")
surv_clim_scalers <- readRDS("../../survival_all_clim_scalers.RDS")
rec_clim_scalers <- readRDS("../../recruitment_all_clim_scalers.RDS")

# Source scripts with functions to import parameters from MCMC iterations
# Returns yearly slopes and intercepts if doYear != NA
source("../vitalRateRegs/survival/import2ipm.R")
source("../vitalRateRegs/growth/import2ipm.R")
source("../vitalRateRegs/recruitment/import2ipm.R")


####
####  LOOP OVER SPECIES AND PROJECT IPM EACH TIME STEP
####
for(spp in 1:length(spp_list)){
  doSpp <- spp_list[spp]
  outfile1 <- paste("./results/forecast_horizon/",doSpp,"_final_year_cover.RDS",sep="")
  sppCode <- which(spp_list==doSpp)
  n_spp <- length(spp_list) # this is needed b/c all 4 spp parameters are imported at once
  
  
  ###
  ### IPM initializations
  ###
  Atotal <- 10000 #area of quadrat in cm^2
  bigM <- c(75,50,10,50)[sppCode] # set matrix dimension for each species; retreive focal spp dim
  maxSize <- c(2500,120,25,100)[sppCode] # set max size for each species; retreive focal spp dim
  L <- log(0.2) # lower limit for size classes
  U <- log(maxSize)*1.1 # upper limit for size classes    
  
  ##  Empty objects for storage -- all set to NULL
  v <- v.r <- b.r <- expv <- Cr <- WmatG <- WmatS <- NULL
  h <- r.L <- r.U <- Ctot <- NULL
  
  # b---boundary points. Note: b chops up the size interval (L-U) into bigM-equal-sized portions.
  b <- L+c(0:bigM)*(U-L)/bigM #log-transformed
  
  # mid points
  v <- 0.5*(b[1:bigM]+b[2:(bigM+1)]) #log-transformed
  
  # step size for midpoint rule. (see equations 4 and 5 in Ellner and Rees (2006) Am Nat.)
  h <- v[2]-v[1] #log-transformed 
  
  # variables for Wr approximation---radius
  b.r <- sqrt(exp(b)/pi) 
  v.r <- sqrt(exp(v)/pi)
  expv <- exp(v) # for size
  r.L <- sqrt(exp(L)/pi) # the lower size limit; radius
  r.U <- sqrt(exp(U)/pi) # the upper size limit; radius
  WmatG <- rep(NA,length(v.r))  # storage of size-specific conspecific W values for each species
  WmatS <- rep(NA,length(v.r))  # storage of size-specific conspecific W values for each species
  
  tmp <- range(v.r) # this is the 'real" range, not log-transformed...
  size.range <- seq(tmp[1],tmp[2],length=bigM) # range across all possible sizes; 'real' radius for each size stage
  
  
  ###
  ### Load Alpha Values for crowding effects (see Chu and Adler 2015)
  ###
  alpha_grow <- read.csv("../../alpha_list_growth.csv")
  alpha_surv <- read.csv("../../alpha_list_survival.csv")
  # Pull out alphas only for the doSpp
  alphaG <- subset(alpha_grow, Site=="Montana")$Alpha[sppCode]
  alphaS <- subset(alpha_surv, Site=="Montana")$Alpha[sppCode]
  
  
  ###
  ### Read in IPM functions
  ###
  source("one_step_ahead_functions.r")
  source("vital_rate_ipm_functions.R")
  
  
  ###
  ### Read in Observation Data
  ###
  # Read in climate data 
  clim_data <- read.csv("../../weather/Climate.csv")
  clim_data <- clim_data[,c("year", "pptLag", "ppt1","ppt2","TmeanSpr1","TmeanSpr2")] # subset and reorder to match regression param import
  
  # Get calendar years
  coverDat <- read.csv("../../speciesData/quadAllCover.csv")
  years <- unique(coverDat$year)
  years <- years[1:(length(years)-1)] #lop off 1945 since no climate for that year
  
  # Get only quads with cover > 0 in year 1944
  tmp<-subset(coverDat,Species==doSpp)
  tmp<-tmp[,c("quad","year","propCover")]
  names(tmp)<-c("quad","t1","obs.cover.t1")
  tmp <- subset(tmp, t1==44)
  quads2sim <- tmp[which(tmp$obs.cover.t1>0),"quad"]
  
  # Read in quad inventory
  quadYearList <- read.csv("../../speciesData/quad_inventory_2.csv")
  quads <- names(quadYearList)
  quads <- quads[which(quads %in% quads2sim)] # skips quads with 0 cover in final year
  
  # Get quad-group info
  quadGroups <- read.csv("../../speciesData/quad_info.csv")
  names(quadGroups) <- c("quad","Group","grazing")
  
  # Read in genet data
  gen_dat <- read.csv(paste("../../speciesData/",doSpp,"/survD.csv",sep=""))
  
  # Read in weather time series
  clim_names <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
  weather <- clim_data[,clim_names]
  weather$inter1 <- weather$ppt1*weather$TmeanSpr1
  weather$inter2 <- weather$ppt2*weather$TmeanSpr2
  
  # Scale climate by means and sds specific to each vital rate
  Gscalers <- subset(growth_clim_scalers, yearout==1 & species==doSpp)
  Sscalers <- subset(surv_clim_scalers, yearout==1 & species==doSpp)
  Rscalers <- subset(rec_clim_scalers, yearout==1 & species==doSpp)
  weatherG <- weatherS <- weatherR <- weather
  for(iclim in 1:nrow(weather)){
    weatherG[iclim,] <- (weatherG[iclim,] - Gscalers[1:length(weatherG[iclim,]),"means"])/Gscalers[1:length(weatherG[iclim,]),"sds"]
    weatherS[iclim,] <- (weatherS[iclim,] - Sscalers[1:length(weatherS[iclim,]),"means"])/Sscalers[1:length(weatherS[iclim,]),"sds"]
    weatherR[iclim,] <- (weatherR[iclim,] - Rscalers[1:length(weatherR[iclim,]),"means"])/Rscalers[1:length(weatherR[iclim,]),"sds"]
  }
  
  
  ###
  ### LOOP THROUGH QUADS AND YEARS ---> PREDICT COVER CHANGE EACH YEAR
  ###
  output_matrix <- as.data.frame(matrix(data=NA, ncol=4,nrow=reps*length(years)*length(quads)))
  counter <- 1 # counter to keep track of looping
  for(iQ in 1:length(quads)){
    # Identify current focal group of quadrats
    doGroup <- quadGroups$Group[which(quadGroups$quad==quads[iQ])]
    doGroup <- as.numeric(doGroup)  # from factor to scalar
    
    for(iYr in 1:length(years)){
      # Identify current year to project
      doYear <- years[iYr]
      
      # Make vector of years to simulate
      years2sim <- (doYear+1):max(years)
      
      # Set recruit size parameters
      rec_size_mean <- numeric(n_spp) # storage vector for mean recruit size
      rec_size_var <- numeric(n_spp) # storage vector for variance in recruit size
      # Loop over species to fill recruit size vectors
      for(i in 1:n_spp){
        infile=paste("../../speciesData/",spp_list[i],"/recSize.csv",sep="")
        recSize=read.csv(infile)
        rec_size_mean[i]=mean(log(recSize$area))
        rec_size_var[i]=var(log(recSize$area))
      }
      
      # Only work with complete transitions
      if(!is.na(quadYearList[iYr,iQ]) & !is.na(quadYearList[iYr+1,iQ])){
        
        tmpD <- subset(gen_dat,quad==quads[iQ] & year==years[iYr])
        nt.init <- table(cut(log(tmpD$area),breaks=b))
        cover.t0 <- sum(tmpD$area)/Atotal
        
        if(cover.t0>0){
          
          # create Cr function for calculating neighborhood crowding
          Ctot <- sum(expv*nt.init) #h*sum(expv*nt.init) #total cover--multiply by h or not?
          Cr <- splinefun(b.r,c(0,cumsum(expv*nt.init)),method="natural") #Cr is a function             
          
          for(iRep in 1:reps){
            nt.new <- nt.init
            for(simyear in years2sim){
              # Get this year's weather
              yearid <- simyear-(min(years)-1) # gets year ID, rather than actual year
              # climyear <- 1900+simyear # tack on 1900 to the doYear to match format in climate data frame
              weather <- list(weatherG[yearid,],weatherS[yearid,],weatherR[yearid,])
              names(weather) <- c("grow_weather", "surv_weather", "rec_weather")
              # Call IPM script ---> No random year effects, climate only
              nt.new <- projectIPM(nt=nt.new,doYear=NA,doGroup,weather,sppCode)
            } # end simulation years loop
            
            cover.t1 <- sum(nt.new*exp(v))/Atotal # convert nt size vector to cover of 1m^2 plot
            # Store predicted final year cover
            output_matrix[counter,] <- c(quads[iQ],doYear,iRep,cover.t1)
            counter <- counter+1
          } # next iRep
        } # end cover.t0>0 if
      } # end is.na(year) if
      
      print(paste(quads[iQ],years[iYr],sep=" "))
      flush.console()
      
    } # next iYr
  } # next iQ 
  
  # Merge with observed data for comparison
  coverDat<-subset(coverDat,Species==doSpp)
  coverDat<-coverDat[,c("quad","year","propCover")]
  names(coverDat)<-c("quad","t1","obs.cover.t1")
  coverDat <- subset(coverDat, t1==44)
  names(output_matrix) <- c("quad","startyear","rep","finalyear_cover")
  output2<-merge(output_matrix,coverDat)
  
  # Get the DF structures correct
  output2$startyear <- as.integer(output2$startyear)
  output2$rep <- as.numeric(output2$rep)
  output2$finalyear_cover <- as.numeric(output2$finalyear_cover)
  
  # Write output
 saveRDS(output2,outfile1)
} # end species loop

