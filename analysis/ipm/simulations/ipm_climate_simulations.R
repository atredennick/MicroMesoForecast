##  Single species IPM with climate covariates
##  This script simulates equilibrium cover



####
####  LOAD LIBRARIES
####
library(boot)
library(mvtnorm)
library(msm)
library(statmod) #"Statistical Modeling" 

####
####  INPUTS
####
doGroup <- NA  # NA for spatial avg., values 1-6 for a specific group
initialCover <- c(0.01) # arbitrarily low starting cover



####
####  SET UP PARAMETERS
####
sppCode <- which(spp_list==doSpp) # gets numeric code for current species (1:4)

##  Fetch recruit size parameters
rec_size_mean <- numeric(n_spp)
rec_size_var <- numeric(n_spp)
for(i in 1:n_spp){
  infile=paste("../../speciesData/",spp_list[i],"/recSize.csv",sep="")
  recSize=read.csv(infile)
  rec_size_mean[i]=mean(log(recSize$area)) # mean new recuit size
  rec_size_var[i]=var(log(recSize$area)) # variance in new recruitm size
}

##  Fetch alpha current species
alphaG <- subset(alpha_grow, Site=="Montana")$Alpha[sppCode]
alphaS <- subset(alpha_surv, Site=="Montana")$Alpha[sppCode]



####
#### SIMULATION LENGTH, MATRIX SIZE, AND INITIAL VECTORS
####
Atotal <- 10000 # Area of 100cm x 100cm quadrat
bigM <- c(75,50,10,50) [sppCode] # Get matrix dimension for current species
maxSize <- c(2500,120,25,100) [sppCode] # Get maximum size for current species    

##  Initialize object names, set to NULL
v <- v.r <- b.r <- expv <- Cr <- WmatG <- WmatS <- NULL 
h <- r.L <- r.U <- Ctot <- NULL

# stuff for numerical approximation....
# minimum (0.2) and maximum sizes (1.1*maximum size from data)
# log-transformation...
# then L, U, b, v, h are all log-transformed....
L <- log(0.2)
U <- log(maxSize)*1.1     
  
# b: boundary points. 
# Note: b chops up the size interval (L-U) into bigM-equal-sized portions.
b <- L+c(0:bigM)*(U-L)/bigM # log-transformed
  
# v: middle points (in Ellner's script, v is y); 
# calculates the middle of each n-equal-sized portion.
v <- 0.5*(b[1:bigM]+b[2:(bigM+1)]) # log-transformed
 
# step size for midpoint rule. 
# see equations 4 and 5 in Ellner and Rees (2006) Am Nat.
h <- v[2]-v[1] # log-transoformed 

# variables for Wr approximation: radius (has to do with crowding)
b.r <- sqrt(exp(b)/pi) 
v.r <- sqrt(exp(v)/pi)
expv <- exp(v) # for size
r.L <- sqrt(exp(L)/pi) # the lower size limit; radius
r.U <- sqrt(exp(U)/pi) # the upper size limit; radius

# storage of size-specific conspecific W values for growth for each species
WmatG <- rep(NA,length(v.r))

# storage of size-specific conspecific W values for survival for each species
WmatS <- rep(NA,length(v.r))

tmp <- range(v.r) # this is the 'real" range, not log-transformed...
size.range <- seq(tmp[1],tmp[2],length=bigM) # range across all possible sizes; 'real' radius for each size stage



####
####  UTILITY FUNCTIONS
####
# Using u and v for size variables instead of x and y, 
# so that x and y can be used for locations in space.
# v above is also for the size (middle point size---area)

#!! Put all component matrices into 3-dimensional arrays !!#

# Combined Kernel
make.K.values <- function(v,u,muWG,muWS, #state variables
                       Rpars,rpa,Gpars,Spars,doYear,doSpp,
                       weatherG,weatherS)  #vital rate arguments
{
  f(v,u,Rpars,rpa,doSpp)+S(u,muWS,Spars,doYear,doSpp,weatherG)*G(v,u,muWG,Gpars,doYear,doSpp,weatherG) 
}

# Function to make iteration matrix based only on mean crowding
make.K.matrix <- function(v,muWG,muWS,Rpars,rpa,Gpars,Spars,doYear,doSpp,weatherG,weatherS) {
  muWS<-rep(muWS,length(v))
  muWG<-rep(muWG,length(v))
  K.matrix <- outer(v,v,make.K.values,muWG,muWS,Rpars,rpa,Gpars,Spars,doYear,doSpp,weatherG,weatherS)
  return(h*K.matrix)
}

# Function to calculate size-dependent crowding, assuming no overlap
wriG <- function(r){
  return(2*pi*integrate(function(z) z*exp(-alphaG*(z^2))*Cr(z-r),r,r+r.U)$value+
         pi*Ctot*exp(-alphaG*((r+r.U)^2))/alphaG)
}
WriG <- Vectorize(wriG,vectorize.args="r") #Wri, and wri

wriS <- function(r){
  return(2*pi*integrate(function(z) z*exp(-alphaS*(z^2))*Cr(z-r),r,r+r.U)$value+
    pi*Ctot*exp(-alphaS*((r+r.U)^2))/alphaS)
}
WriS <- Vectorize(wriS,vectorize.args="r") #Wri, and wri

# Function to sum total cover of each species
sumCover <- function(v,nt,h,Atotal){ 
  out <- h*sum(nt*exp(v))/Atotal 
  return (out)
}

# Function to sum total density of each species
sumN <- function(nt,h){ 
  out <- h*sum(nt)
  return (out)
} 

# Function to calculate initial nt vector, given initial cover
initial.nt <- function(initialCover,Atotal,h,v){ # here v is a list with one element
  tmp <- 0
  for(i in 1:length(v)){
    tmp <- tmp+exp(v[i])
    tmp <- tmp
  }
  nt <- initialCover*Atotal/(h*tmp)
  return (nt)
}



####
####  RUN EQUILIBRIUM SIMULATION
####
# Initial population density vector for each size 
# (so need integral when for total cover & density)
initial_nt <- initial.nt(initialCover,Atotal,h,v)
nt <- matrix(initial_nt,bigM) # to store the population size
new.nt <- nt

# ntTicks is for recording the size vector for every simulation time
sizeSave <- matrix(NA,bigM,tlimit)
sizeSave[,1] <- nt

# set up vector to record cover
covSave <- rep(NA,tlimit) 
covSave[1] <- sumCover(v,nt,h,Atotal)

# set up vector to record density 
Nsave <- rep(NA,tlimit) 
Nsave[1] <- sumN(nt,h)

for (t in 2:(tlimit)){
  doYear <- NA # no random year effects; climate only
  doClim <- climYr[t]-(min(climYr)-1)
  weather <- clim_data[clim_data$year==(1900+climYr[t]),2:6]
  weather$inter1 <- weather$ppt1*weather$TmeanSpr1
  weather$inter2 <- weather$ppt2*weather$TmeanSpr2
  
  # Scale climate by means and sds specific to each vital rate
  weatherG <- (weather - Gscalers[1:length(weather),"means"])/Gscalers[1:length(weather),"sds"]
  weatherS <- (weather - Sscalers[1:length(weather),"means"])/Sscalers[1:length(weather),"sds"]
  weatherR <- (weather - Rscalers[1:length(weather),"means"])/Rscalers[1:length(weather),"sds"]
  
  # get vital rate parameters
  Spars<-getSurvCoefs(doYear,doGroup)
  Gpars<-getGrowCoefs(doYear,doGroup)
  Rpars<-getRecCoefs(doYear,doGroup)
  Rpars$sizeMean <- rec_size_mean
  Rpars$sizeVar <- rec_size_var
  
  # calculate intraspecific crowding 
  Ctot=h*sum(expv*nt) 
  Cr=splinefun(b.r,h*c(0,cumsum(expv*nt)),method="natural") #Cr is a function      
  WfunG=splinefun(size.range,WriG(size.range))
  WmatG=WfunG(v.r)/Atotal
  WfunS=splinefun(size.range,WriS(size.range))
  WmatS=WfunS(v.r)/Atotal
  
  # get recruits per area
  cover<-covSave[t-1]
  rpa<-rep(0,n_spp)
  rpa[sppCode]=get_rpa(Rpars,cover,weatherR)  # set up this way to eventually run multiple species at once
  
  # make kernel and project population
  K.matrix=make.K.matrix(v,WmatG,WmatS,Rpars,rpa,Gpars,Spars,doYear,sppCode,weatherG,weatherS)  
  new.nt=K.matrix%*%nt
  sizeSave[,t]=new.nt/sum(new.nt)  
  
  # update and store state variables
  nt=new.nt
  covSave[t]<-sumCover(v,nt,h,Atotal)  # store the cover as cm^2/cm^2
  Nsave[t]<-sumN(nt,h)
  
#   if(t > 1 & covSave[t]==0) {
#     initial_nt<-initial.nt(covSave[t-2],Atotal,h,v)
#     nt<-matrix(initial_nt,bigM)
#     covSave[i] <- sumCover(v,nt,h,Atotal)
#   }
  if(sum(is.na(nt))>0) browser() # check for errors

  print(t) # print the simulation steps
  flush.console()
} # next time step

