## This script loops through all observed year-to-year transitions,
## sets up initial conditions, and then calls the IPM to project the population one time step
##
##  Author: Peter Adler and Andrew Tredennick
##  Updated: 4-7-2016

##  Clear the Workspace
rm(list=ls())

##  Set Working Directory
root <- ifelse(.Platform$OS.type=="windows","c:/repos","~/repos"); # modify as needed
setwd(paste(root,"/MicroMesoForecast/montana/ipm/simulations",sep="")); # modify as needed 



####
####  PRELIMINARIES
####
spp_list <- c("BOGR","HECO","PASM","POSE") # all Montana species
reps <- 2  # number of times to simulate each quadrat x year transition



####
####  LOOP OVER SPECIES AND PROJECT IPM EACH TIME STEP
####
for(spp in 1:length(spp_list)){
  doSpp <- spp_list[spp]
  outfile1 <- paste("./one_step_results/",doSpp,"_sim_cover_1step_ahead.csv",sep="")
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
  
  # Read in quad inventory
  quadYearList <- read.csv("../../speciesData/quad_inventory_2.csv")
  quads <- names(quadYearList)
  
  # Get quad-group info
  quadGroups <- read.csv("../../speciesData/quad_info.csv")
  names(quadGroups) <- c("quad","Group","grazing")
  
  # Read in genet data
  gen_dat <- read.csv(paste("../../speciesData/",doSpp,"/survD.csv",sep=""))
  
  
  ###
  ### LOOP THROUGH QUADS AND YEARS ---> PREDICT COVER CHANGE EACH YEAR
  ###
  output <- data.frame(quad="c",t0=1,t1=2,rep=1,
                       cover.t0=1,cover.t1=1,
                       stringsAsFactors=F) # storage data frame for predictions
  counter <- 0 # counter to keep track of looping
  for(iQ in 1:length(quads)){
    # Identify current focal group of quadrats
    doGroup <- quadGroups$Group[which(quadGroups$quad==quads[iQ])]
    doGroup <- as.numeric(doGroup)  # from factor to scalar
    
    for(iYr in 1:length(years)){
      # Identify current year to project
      doYear <- years[iYr]
      
      # Source scripts with functions to import parameters from MCMC iterations
      # Returns yearly slopes and intercepts if doYear != NA
      source("../vitalRateRegs/validation/survival/import2ipm.R")
      source("../vitalRateRegs/validation/growth/import2ipm.R")
      source("../vitalRateRegs/validation/recruitment/import2ipm.R")
      
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
      
      # Get this year's weather
      yearid <- doYear-(min(years)-1) # gets year ID, rather than actual year
      climyear <- 1900+doYear # tack on 1900 to the doYear to match format in climate data frame
      weather <- clim_data[clim_data$year==climyear,2:6]
      weather$inter1 <- weather$ppt1*weather$TmeanSpr1
      weather$inter2 <- weather$ppt2*weather$TmeanSpr2
      
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
            # Call IPM script ---> No random year effects, climate only
            nt.new <- projectIPM(nt=nt.init,doYear=NA,doGroup,weather,sppCode)
            cover.t1 <- sum(nt.new*exp(v))/Atotal # convert nt size vector to cover of 1m^2 plot
            
              # Store predicted cover
            counter <- counter+1
            output[counter,] <- NA
            output$quad[counter] <- quads[iQ]
            output[counter,2:6] <- c(doYear,doYear+1,iRep,cover.t0,cover.t1)
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
  output2<-merge(output,coverDat)
  
  # Write output
  write.csv(output2,outfile1,row.names=F)
} # end species loop




