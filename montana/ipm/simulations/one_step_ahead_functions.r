
#--------------------------------------
# VITAL RATE FUNCTIONS AND PARAMETERS
#--------------------------------------

source("vital_rate_ipm_functions.R")

# get calendar years
raw_data <- read.csv("../../speciesData/quadAllCover.csv")
years <- unique(raw_data$year)
years <- years[1:(length(years)-1)] #lop off 1945 since no climate for that year

source("../vitalRateRegs/survival/import2ipm.R")
source("../vitalRateRegs/growth/import2ipm.R")
source("../vitalRateRegs/recruitment/import2ipm.R")

# get recruit size parameters
rec_size_mean <- numeric(n_spp)
rec_size_var <- numeric(n_spp)
for(i in 1:n_spp){
  infile=paste("../../speciesData/",spp_list[i],"/recSize.csv",sep="")
  recSize=read.csv(infile)
  rec_size_mean[i]=mean(log(recSize$area))
  rec_size_var[i]=var(log(recSize$area))
}

# get alphas values (needed to calculate neighborhood crowding)
alpha_grow <- read.csv("../../alpha_list_growth.csv")
alpha_surv <- read.csv("../../alpha_list_survival.csv")
# pull out alphas only for the doSpp
alphaG <- subset(alpha_grow, Site=="Montana")$Alpha[sppCode]
alphaS <- subset(alpha_surv, Site=="Montana")$Alpha[sppCode]

#============================================================================================#
# (III) Utility functions
#============================================================================================#

#Using u and v for size variables instead of x and y, so that x and y can be used for locations in space
#v above is also for the size (middle point size---area)
# load the necessary libraries
library(boot)
library(mvtnorm)
library(msm)
library(statmod) #"Statistical Modeling" 

##### Put all component matrices into 3-dimensional arrays 

# combined kernel
make.K.values=function(v,u,muWG,muWS, #state variables
                       Rpars,rpa,Gpars,Spars,doYear,doSpp,weather)  #vital rate arguments
{
  f(v,u,Rpars,rpa,doSpp)+S(u,muWS,Spars,doYear,doSpp,weather)*G(v,u,muWG,Gpars,doYear,doSpp,weather) 
}

# Function to make iteration matrix based only on mean crowding
make.K.matrix=function(v,muWG,muWS,Rpars,rpa,Gpars,Spars,doYear,doSpp,weather) {
  #muWS=expandW(v,v,muWS)  # for multispecies models
  #muWG=expandW(v,v,muWG)
  muWS<-rep(muWS,length(v))
  muWG<-rep(muWG,length(v))
  K.matrix=outer(v,v,make.K.values,muWG,muWS,Rpars,rpa,Gpars,Spars,doYear,doSpp,weather)
  return(h*K.matrix)
}

# next function only needed for multispecies models
# # Function to format the W matrix for the outer product
# expandW=function(v,u,W){
#    if(dim(W)[1]!=length(u)) stop("Check size of W")
#    Nspp=dim(W)[2]
#    W=as.vector(W)
#    W=matrix(W,length(W),ncol=length(v))
#    W=as.vector(t(W))
#    W=matrix(W,nrow=length(u)*length(v),ncol=Nspp)
#    return(W)
# }

# Function to calculate size-dependent crowding, assuming no overlap
wriG=function(r){
  return(2*pi*integrate(function(z) z*exp(-alphaG*(z^2))*Cr(z-r),r,r+r.U)$value+
         pi*Ctot*exp(-alphaG*((r+r.U)^2))/alphaG)
}
WriG=Vectorize(wriG,vectorize.args="r") #Wri, and wri

wriS=function(r){
  return(2*pi*integrate(function(z) z*exp(-alphaS*(z^2))*Cr(z-r),r,r+r.U)$value+
    pi*Ctot*exp(-alphaS*((r+r.U)^2))/alphaS)
}
WriS=Vectorize(wriS,vectorize.args="r") #Wri, and wri

# Function to sum total cover of each species
sumCover<-function(v,nt,h,Atotal){ 
  out<-h*sum(nt*exp(v))/Atotal 
  return (out)
}

# Function to sum total density of each species
sumN<-function(nt,h){ 
  out<-h*sum(nt)
  return (out)
} 
     
# Function to do an image plot of a matrix in the usual orientation, A(1,1) at top left  
matrix.image<-function(x,y,A,col=topo.colors(100),...) {
  nx<-length(x)
  ny<-length(y) 
	x1<-c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]) 
	y1<-c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]) 
	image(list(x=x,y=y,z=t(A)),xlim=x1,ylim=rev(y1),col=col,bty="u",...)  
}
