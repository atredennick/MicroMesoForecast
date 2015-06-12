
# single species IPM with climate covariates
# this script simulates equilibrium cover

#============================================================
# (I) INPUTS
#============================================================
doGroup=NA  # NA for spatial avg., values 1-6 for a specific group
initialCover<-c(0.01)
# nMCMC<-3000 # max number of MCMC iterations to draw parameters from

#============================================================
# (II) LOAD VITAL RATE FUNCTIONS & SET UP PARAMETERS
#============================================================

sppCode<-which(spp_list==doSpp)
# n_spp<-length(spp_list) # this is needed b/c all 4 spp parameters are imported at once

# get recruit size parameters
rec_size_mean <- numeric(n_spp)
rec_size_var <- numeric(n_spp)
for(i in 1:n_spp){
  infile=paste("../../speciesData/",spp_list[i],"/recSize.csv",sep="")
  recSize=read.csv(infile)
  rec_size_mean[i]=mean(log(recSize$area))
  rec_size_var[i]=var(log(recSize$area))
}

# pull out alphas only for the doSpp
alphaG <- subset(alpha_grow, Site=="Montana")$Alpha[sppCode]
alphaS <- subset(alpha_surv, Site=="Montana")$Alpha[sppCode]

#============================================================================================#
# (II) Simulation length, Matrix size and initial vectors
#============================================================================================#

Atotal=10000 #Area of 100cm x 100cm quadrat

bigM=c(75,50,10,50) [sppCode]            #Set matrix dimension for each species
maxSize=c(2500,120,25,100) [sppCode]    

v=v.r=b.r=expv=Cr=WmatG=WmatS=NULL #a LIST for different species
h=r.L=r.U=Ctot=NULL

# stuff for numerical approximation....
  # minimum (0.9*minimum size from data) and maximum sizes (1.1*maximum size from data)
  # log-transformation...
  # then L, U, b, v, h are all log-transformed....
L=log(0.2)
U=log(maxSize)*1.1     
  
# b---boundary points. Note: b chops up the size interval (L-U) into bigM-equal-sized portions.
b = L+c(0:bigM)*(U-L)/bigM #log-transformed
  
# v---middle points (in Ellner's script, v is y); calculates the middle of each n-equal-sized portion.
v = 0.5*(b[1:bigM]+b[2:(bigM+1)]) #log-transformed
 
# step size for midpoint rule. (see equations 4 and 5 in Ellner and Rees (2006) Am Nat.)
h = v[2]-v[1] #log-transoformed 

# variables for Wr approximation---radius
b.r=sqrt(exp(b)/pi) 
v.r=sqrt(exp(v)/pi)
expv=exp(v) #for size
r.L=sqrt(exp(L)/pi) #the lower size limit; radius
r.U=sqrt(exp(U)/pi) #the upper size limit; radius
WmatG=rep(NA,length(v.r))  # storage of size-specific conspecific W values for each species
WmatS=rep(NA,length(v.r))  # storage of size-specific conspecific W values for each species

tmp=range(v.r)#this is the 'real" range, not log-transformed...
size.range=seq(tmp[1],tmp[2],length=bigM) # range across all possible sizes; 'real' radius for each size stage

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
                       Rpars,rpa,Gpars,Spars,doYear,doSpp,weather_surv, weather_grow)  #vital rate arguments
{
  f(v,u,Rpars,rpa,doSpp)+S(u,muWS,Spars,doYear,doSpp,weather_surv)*G(v,u,muWG,Gpars,doYear,doSpp,weather_grow) 
}

# Function to make iteration matrix based only on mean crowding
make.K.matrix=function(v,muWG,muWS,Rpars,rpa,Gpars,Spars,doYear,doSpp,weather_surv, weather_grow) {
  #muWS=expandW(v,v,muWS)  # for multispecies models
  #muWG=expandW(v,v,muWG)
  muWS<-rep(muWS,length(v))
  muWG<-rep(muWG,length(v))
  K.matrix=outer(v,v,make.K.values,muWG,muWS,Rpars,rpa,Gpars,Spars,doYear,doSpp,weather_surv, weather_grow)
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

###  calculate the initial nt value, given an inital cover 
initial.nt<-function(initialCover,Atotal,h,v){ ## here v is a list with one element
  tmp<-0
  for(i in 1:length(v)){
    tmp<-tmp+exp(v[i])
    tmp<-tmp
  }
  nt<-initialCover*Atotal/(h*tmp)
  return (nt)
}

#============================================================================================#
# (IV) Calculate the equilibrium areas.
#============================================================================================# 

# initial population density vector for each size (so need integral when for total cover & density)
initial_nt<-initial.nt(initialCover,Atotal,h,v)
nt<-matrix(initial_nt,bigM) #to store the population size
new.nt<-nt

# ntTicks is for recording the size vector for every simulation time
sizeSave<-matrix(NA,bigM,tlimit)
sizeSave[,1]<-nt

# set up vector to record cover
covSave = rep(NA,tlimit) 
covSave[1] <- sumCover(v,nt,h,Atotal)

# set up vector to record density 
Nsave <- rep(NA,tlimit) 
Nsave[1]<-sumN(nt,h)

for (i in 2:(tlimit)){
  doYear<-yrSave[i]-(min(yrSave)-1)
  doClim <- climYr[i]-(min(climYr)-1)
  weather<-clim_data[clim_data$year==(1900+climYr[i]),2:6]
  weather$inter1 <- weather$ppt1*weather$TmeanSpr1
  weather$inter2 <- weather$ppt2*weather$TmeanSpr2
  
  weather_grow <- weather
  weather_surv <- weather
  weather_rec <- weather
  
  #Reset weather for vital rate sensitivity
  if(doPpt == "none"){
    weather_grow <- weather_grow
    weather_surv <- weather_surv
    weather_rec <- weather_rec
  }
  if(doTemp == "none"){
    weather_grow <- weather_grow
    weather_surv <- weather_surv
    weather_rec <- weather_rec
  }
  if(doBoth == "none"){
    weather_grow <- weather_grow
    weather_surv <- weather_surv
    weather_rec <- weather_rec
  }
  if(doPpt == "growth"){
    weather_grow <- clim_ppt[clim_ppt$year==(1900+climYr[i]),2:6]
    weather_grow$inter1 <- weather_grow$ppt1*weather_grow$TmeanSpr1
    weather_grow$inter2 <- weather_grow$ppt2*weather_grow$TmeanSpr2
  }
  if(doPpt == "survival"){
    weather_surv <- clim_ppt[clim_ppt$year==(1900+climYr[i]),2:6]
    weather_surv$inter1 <- weather_surv$ppt1*weather_surv$TmeanSpr1
    weather_surv$inter2 <- weather_surv$ppt2*weather_surv$TmeanSpr2
  }
  if(doPpt == "recruitment"){
    weather_rec <- clim_ppt[clim_ppt$year==(1900+climYr[i]),2:6]
    weather_rec$inter1 <- weather_rec$ppt1*weather_rec$TmeanSpr1
    weather_rec$inter2 <- weather_rec$ppt2*weather_rec$TmeanSpr2
  }
  
  if(doPpt == "growth_surv"){
    weather_grow <- clim_ppt[clim_ppt$year==(1900+climYr[i]),2:6]
    weather_grow$inter1 <- weather_grow$ppt1*weather_grow$TmeanSpr1
    weather_grow$inter2 <- weather_grow$ppt2*weather_grow$TmeanSpr2
    weather_surv <- clim_ppt[clim_ppt$year==(1900+climYr[i]),2:6]
    weather_surv$inter1 <- weather_surv$ppt1*weather_surv$TmeanSpr1
    weather_surv$inter2 <- weather_surv$ppt2*weather_surv$TmeanSpr2
  }
  if(doPpt == "growth_rec"){
    weather_grow <- clim_ppt[clim_ppt$year==(1900+climYr[i]),2:6]
    weather_grow$inter1 <- weather_grow$ppt1*weather_grow$TmeanSpr1
    weather_grow$inter2 <- weather_grow$ppt2*weather_grow$TmeanSpr2
    weather_rec <- clim_ppt[clim_ppt$year==(1900+climYr[i]),2:6]
    weather_rec$inter1 <- weather_rec$ppt1*weather_rec$TmeanSpr1
    weather_rec$inter2 <- weather_rec$ppt2*weather_rec$TmeanSpr2
  }
  if(doPpt == "surv_rec"){
    weather_surv <- clim_ppt[clim_ppt$year==(1900+climYr[i]),2:6]
    weather_surv$inter1 <- weather_surv$ppt1*weather_surv$TmeanSpr1
    weather_surv$inter2 <- weather_surv$ppt2*weather_surv$TmeanSpr2
    weather_rec <- clim_ppt[clim_ppt$year==(1900+climYr[i]),2:6]
    weather_rec$inter1 <- weather_rec$ppt1*weather_rec$TmeanSpr1
    weather_rec$inter2 <- weather_rec$ppt2*weather_rec$TmeanSpr2
  }
  if(doTemp == "growth"){
    weather_grow <- clim_temp[clim_temp$year==(1900+climYr[i]),2:6]
    weather_grow$inter1 <- weather_grow$ppt1*weather_grow$TmeanSpr1
    weather_grow$inter2 <- weather_grow$ppt2*weather_grow$TmeanSpr2
  }
  if(doTemp == "survival"){
    weather_surv <- clim_temp[clim_temp$year==(1900+climYr[i]),2:6]
    weather_surv$inter1 <- weather_surv$ppt1*weather_surv$TmeanSpr1
    weather_surv$inter2 <- weather_surv$ppt2*weather_surv$TmeanSpr2
  }
  if(doTemp == "recruitment"){
    weather_rec <- clim_temp[clim_temp$year==(1900+climYr[i]),2:6]
    weather_rec$inter1 <- weather_rec$ppt1*weather_rec$TmeanSpr1
    weather_rec$inter2 <- weather_rec$ppt2*weather_rec$TmeanSpr2
  }
  if(doTemp == "growth_surv"){
    weather_grow <- clim_temp[clim_temp$year==(1900+climYr[i]),2:6]
    weather_grow$inter1 <- weather_grow$ppt1*weather_grow$TmeanSpr1
    weather_grow$inter2 <- weather_grow$ppt2*weather_grow$TmeanSpr2
    weather_surv <- clim_temp[clim_temp$year==(1900+climYr[i]),2:6]
    weather_surv$inter1 <- weather_surv$ppt1*weather_surv$TmeanSpr1
    weather_surv$inter2 <- weather_surv$ppt2*weather_surv$TmeanSpr2
  }
  if(doTemp == "growth_rec"){
    weather_grow <- clim_temp[clim_temp$year==(1900+climYr[i]),2:6]
    weather_grow$inter1 <- weather_grow$ppt1*weather_grow$TmeanSpr1
    weather_grow$inter2 <- weather_grow$ppt2*weather_grow$TmeanSpr2
    weather_rec <- clim_temp[clim_temp$year==(1900+climYr[i]),2:6]
    weather_rec$inter1 <- weather_rec$ppt1*weather_rec$TmeanSpr1
    weather_rec$inter2 <- weather_rec$ppt2*weather_rec$TmeanSpr2
  }
  if(doTemp == "surv_rec"){
    weather_surv <- clim_temp[clim_temp$year==(1900+climYr[i]),2:6]
    weather_surv$inter1 <- weather_surv$ppt1*weather_surv$TmeanSpr1
    weather_surv$inter2 <- weather_surv$ppt2*weather_surv$TmeanSpr2
    weather_rec <- clim_temp[clim_temp$year==(1900+climYr[i]),2:6]
    weather_rec$inter1 <- weather_rec$ppt1*weather_rec$TmeanSpr1
    weather_rec$inter2 <- weather_rec$ppt2*weather_rec$TmeanSpr2
  }
  if(doBoth=="growth"){
    weather_grow <- clim_both[climt_both$year==(1900+climYr[i]),2:6]
    weather_grow$inter1 <- weather_grow$ppt1*weather_grow$TmeanSpr1
    weather_grow$inter2 <- weather_grow$ppt2*weather_grow$TmeanSpr2
  }
  if(doBoth=="survival"){
    weather_surv <- clim_both[clim_both$year==(1900+climYr[i]),2:6]
    weather_surv$inter1 <- weather_surv$ppt1*weather_surv$TmeanSpr1
    weather_surv$inter2 <- weather_surv$ppt2*weather_surv$TmeanSpr2
  }
  if(doBoth=="recruitment"){
    weather_rec <- clim_both[clim_both$year==(1900+climYr[i]),2:6]
    weather_rec$inter1 <- weather_rec$ppt1*weather_rec$TmeanSpr1
    weather_rec$inter2 <- weather_rec$ppt2*weather_rec$TmeanSpr2
  }
  if(doBoth=="growth_surv"){
    weather_grow <- clim_both[climt_both$year==(1900+climYr[i]),2:6]
    weather_grow$inter1 <- weather_grow$ppt1*weather_grow$TmeanSpr1
    weather_grow$inter2 <- weather_grow$ppt2*weather_grow$TmeanSpr2
    weather_surv <- clim_both[clim_both$year==(1900+climYr[i]),2:6]
    weather_surv$inter1 <- weather_surv$ppt1*weather_surv$TmeanSpr1
    weather_surv$inter2 <- weather_surv$ppt2*weather_surv$TmeanSpr2
  }
  if(doBoth=="growth_rec"){
    weather_grow <- clim_both[climt_both$year==(1900+climYr[i]),2:6]
    weather_grow$inter1 <- weather_grow$ppt1*weather_grow$TmeanSpr1
    weather_grow$inter2 <- weather_grow$ppt2*weather_grow$TmeanSpr2
    weather_rec <- clim_both[clim_both$year==(1900+climYr[i]),2:6]
    weather_rec$inter1 <- weather_rec$ppt1*weather_rec$TmeanSpr1
    weather_rec$inter2 <- weather_rec$ppt2*weather_rec$TmeanSpr2
  }
  if(doBoth=="surv_rec"){
    weather_surv <- clim_both[clim_both$year==(1900+climYr[i]),2:6]
    weather_surv$inter1 <- weather_surv$ppt1*weather_surv$TmeanSpr1
    weather_surv$inter2 <- weather_surv$ppt2*weather_surv$TmeanSpr2
    weather_rec <- clim_both[clim_both$year==(1900+climYr[i]),2:6]
    weather_rec$inter1 <- weather_rec$ppt1*weather_rec$TmeanSpr1
    weather_rec$inter2 <- weather_rec$ppt2*weather_rec$TmeanSpr2
  }
  
  # get vital rate parameters
  doYear=NA #no random year effects, just climate variation
  Spars<-getSurvCoefs(doYear,doGroup)
  Gpars<-getGrowCoefs(doYear,doGroup)
  Rpars<-getRecCoefs(doYear,doGroup)
  Rpars$sizeMean <- rec_size_mean
  Rpars$sizeVar <- rec_size_var
  
#   Gpars$clim[] <- 0
#   Gpars$slopeXclim[] <- 0
#   Spars$clim[] <- 0
#   Spars$slopeXclim[] <- 0
#   Rpars$clim[] <- 0
  
  # calculate intraspecific crowding 
  Ctot=h*sum(expv*nt) 
  Cr=splinefun(b.r,h*c(0,cumsum(expv*nt)),method="natural") #Cr is a function      
  WfunG=splinefun(size.range,WriG(size.range))
  WmatG=WfunG(v.r)/Atotal
  WfunS=splinefun(size.range,WriS(size.range))
  WmatS=WfunS(v.r)/Atotal
  
  # get recruits per area
  cover<-covSave[i-1]
  rpa<-rep(0,n_spp); rpa[sppCode]=get_rpa(Rpars,cover,weather_rec)  # set up this way to eventually run multiple species at once
  
  # make kernel and project population
  K.matrix=make.K.matrix(v,WmatG,WmatS,Rpars,rpa,Gpars,Spars,doYear,sppCode,weather_surv, weather_grow)  
  new.nt=K.matrix%*%nt
  sizeSave[,i]=new.nt/sum(new.nt)  
  
  # update and store state variables
  nt=new.nt
  covSave[i]<-sumCover(v,nt,h,Atotal)  # store the cover as cm^2/cm^2
  Nsave[i]<-sumN(nt,h)
  
  if(sum(is.na(nt))>0) browser() # check for errors

  print(i)# print the simulation steps
  flush.console()
  
} # next time step


#============================================================================================#
# (V) Output
#============================================================================================# 
# covSave <- covSave[which(covSave<1)] #this is just for graphical purposes
## Figures
# par(mfrow=c(2,2),tcl=-0.2,mgp=c(2,0.5,0)) 
# plot(covSave*100, type="l")
# plot(burn.in:tlimit,100*covSave[burn.in:tlimit],type="l",xlab="Time",ylab="Cover (%)")
# plot(burn.in:tlimit,Nsave[burn.in:tlimit],type="l",xlab="Time",ylab="Density") 
# plot(1,1,type="n",xlim=c(log(0.15),log(max(maxSize))),ylim=c(0,0.1),xlab="Size",ylab="Frequency")
# lines(v,rowMeans(sizeSave[,(burn.in+1):tlimit])) # average size distribution 
# plot(density(100*covSave[burn.in:tlimit],na.rm = TRUE), xlim=c(0,100))

## Write data tables
output1<-data.frame("time"=burn.in:tlimit,"cover"=covSave[burn.in:tlimit])
output1$species <- sppCode
output1$climsim <- simcode
output1$vital <- vitalcode
output2<-data.frame("time"=burn.in:tlimit,"density"=Nsave[burn.in:tlimit])
output3<-data.frame("size"=v,"frequency"=rowMeans(sizeSave[,(burn.in+1):tlimit]))

write.table(output1,outfile1,row.names=F,sep=",")
write.table(output2,outfile2,row.names=F,sep=",")
write.table(output3,outfile3,row.names=F,sep=",")

