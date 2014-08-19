# Multispecies, spatially implicit IPM
# This version makes it possible to assume "no overlap"
# for intraspecific competition only or intra- and interspecific competition

# PBA 8/17/10

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

outfile1="ipm_cover_T10.csv"
outfile2="stable_size.csv"
obsClimateFile="Climate.csv"
perturbPpt=T
perturbTemp=T
climYrSave=read.csv("climYears.csv")  # use same sequence of years used for observed run
randYrSave=read.csv("randYears.csv")
A=10000 #Area of 100cm x 100cm quadrat
tlimit=2250 ## number of years to simulate
burn.in=250    # years to cut before calculations
sppList=c("ARTR","HECO","POSE","PSSP")
bigM=c(75,75,50,50)     #Set matrix dimension for each species
maxSize=c(3000,202,260,225)    # in cm^2: PSSP=225 HECO=202  POSE=260  ARTR=3000  # minSize=0.2  cm^2
Nyrs=22
doGroup=NA  # NA for spatial avg., values 1-6 for a specific group
constant=F                

#============================================================
# (I) LOAD VITAL RATE PARAMETERS & FUNCTIONS
#============================================================
Nspp=length(sppList)

# get climate data
obsD=read.csv(obsClimateFile)
climD=obsD[,2:NCOL(obsD)]  #drop year column
climD$ppt1.TmeanSpr1=climD$ppt1*climD$TmeanSpr1
climD$ppt2.TmeanSpr2=climD$ppt2*climD$TmeanSpr2
climDmean=colMeans(climD); climDsd=apply(climD,MARGIN=2,FUN=sd)

# set up survival parameters and function
source("survival/import2ipm_noOverlap.r", local=TRUE)
# set up growth parameters and function
source("growth/import2ipm_noOverlap.r", local=TRUE)
# set up recruitment parameters and function
source("recruitment/import2ipm.r", local=TRUE)


# model spatial group variation (or not)
if(!is.na(doGroup)){
  Spars$intcpt=Spars$intcpt+Spars$intcpt.gr[doGroup,]
  Gpars$intcpt=Gpars$intcpt+Gpars$intcpt.gr[doGroup,]
  Rpars$intcpt.yr=Rpars$intcpt.yr+matrix(Rpars$intcpt.gr[doGroup,],Nyrs,Nspp,byrow=T)
}

# PERTURB PARAMETERS -------------------------------------

# perturb precip (increase 1 %)
if(perturbPpt==T){
  pptVars=grep("ppt",names(climD))
  tmp1=0.01*colMeans(climD)
  tmp1=matrix(tmp1,NROW(climD),length(tmp1),byrow=T)
  climD[,pptVars]=climD[,pptVars]+tmp1[,pptVars]
  
  #downscaled climate projections
  #climD$pptLag=climD$pptLag*1.005
  #climD$ppt1=climD$ppt1*1.027
  #climD$ppt2=climD$ppt2*1.027
}

# perturb temperature (increase all temps by 1 % of spring temps)
if(perturbTemp==T){  
  tempVars=grep("Tmean",names(climD))
  tmp1=0.01*mean(climD$TmeanSpr1)
  #tmp1=0.254*mean(climD$TmeanSpr1)   # downscaled projections
  tmp1=matrix(tmp1,NROW(climD),NCOL(climD),byrow=T)
  climD[,tempVars]=climD[,tempVars]+tmp1[,tempVars]
  
}
# recalculate climate interactions
climD$ppt1.TmeanSpr1=climD$ppt1*climD$TmeanSpr1
climD$ppt2.TmeanSpr2=climD$ppt2*climD$TmeanSpr2


if(constant==T){
  # turn off climate variation  
    climD=matrix(climDmean,Nyrs,NCOL(climD),byrow=T) 
  #turn off random year effects  
    #Rpars$intcpt.yr=matrix(Rpars$intcpt.mu,Nyrs,Nspp,byrow=T)
    #Gpars$intcpt.yr[]=0;Gpars$slope.yr[]=0
    #Spars$intcpt.yr[]=0;Spars$slope.yr[]=0    
}

#============================================================================================#
# (II) Simulation length, Matrix size and initial vectors
#============================================================================================#

v=v.r=b.r=expv=Cr=Wmat=list(4)
h=r.L=r.U=Ctot=numeric(4)
for(i in 1:Nspp){
  
  # minimum (0.9*minimum size from data) and maximum sizes (1.1*maximum size from data)
  L=log(0.2)
  U=log(maxSize[i])*1.1     
  
  # boundary points b and mesh points y. Note: b chops up the size interval (L-U) into bigM-equal-sized portions.
  b = L+c(0:bigM[i])*(U-L)/bigM[i] 
  
  # v calculates the middle of each n-equal-sized portion.
  v[[i]] = 0.5*(b[1:bigM[i]]+b[2:(bigM[i]+1)])
  
  # step size for midpoint rule. (see equations 4 and 5 in Ellner and Rees (2006) Am Nat.)
  h[i] = v[[i]][2]-v[[i]][1]  
  
  # variables for Wr approximation
  b.r[[i]]=sqrt(exp(b)/pi)
  v.r[[i]]=sqrt(exp(v[[i]])/pi)
  expv[[i]]=exp(v[[i]])
  r.L[i] = sqrt(exp(L)/pi); 
  r.U[i] = sqrt(exp(U)/pi); 
  Wmat[[i]]=matrix(NA,length(v.r[[i]]),Nspp)  # storage of size-specific W values for each focal species
  
} # next species
tmp=range(v.r)
size.range=seq(tmp[1],tmp[2],length=50) # range across all possible sizes

#============================================================================================#
# (III) Utility functions
#============================================================================================#

# load the necessary libraries
library(boot)
library(mvtnorm)
library(msm)
library(statmod)  

## combined kernel
make.K.values=function(v,u,muW, #state variables
  Rpars,rpa,Gpars,Spars,weather,doYear,doSpp)  #growth arguments
  {
  f(v,u,Rpars,rpa,doSpp)+S(u,muW,Spars,weather,doYear,doSpp)*G(v,u,muW,Gpars,weather,doYear,doSpp) 
}

# Function to make iteration matrix based only on mean crowding
make.K.matrix=function(v,muW,Rpars,rpa,Gpars,Spars,weather,doYear,doSpp) {
       muW=expandW(v,v,muW)
       K.matrix=outer(v,v,make.K.values,muW,Rpars,rpa,Gpars,Spars,weather,doYear,doSpp)
       return(h[doSpp]*K.matrix)
}

# Function to format the W matrix for the outer product
expandW=function(v,u,W){
   if(dim(W)[1]!=length(u)) stop("Check size of W")
   Nspp=dim(W)[2]
   W=as.vector(W)
   W=matrix(W,length(W),ncol=length(v))
   W=as.vector(t(W))
   W=matrix(W,nrow=length(u)*length(v),ncol=Nspp)
   return(W)
}

# Function to calculate size-dependent crowding, assuming no overlap
wrij=function(r,i,j) {
   return(2*pi*integrate(function(z) z*exp(-alpha[i,j]*(z^2))*Cr[[j]](z-r), r, r + r.U[j])$value
     + pi*Ctot[j]*exp(-alpha[i,j]*((r+r.U[j])^2))/alpha[i,j]); 	
}
Wrij=Vectorize(wrij,vectorize.args="r")


# Function to sum total cover of each species
sumCover=function(v,nt,h,A){
   out=lapply(1:Nspp,function(i,v,nt,h,A) h[i]*sum(nt[[i]]*exp(v[[i]]))/A,v=v,nt=nt,h=h,A=A)
   return(unlist(out))
} 

# Function to sum total density of each species
sumN=function(nt,h){
   out=lapply(1:Nspp,function(i,nt,h) h[i]*sum(nt[[i]]),nt=nt,h=h)
   return(unlist(out))
}

# Function to calculate size variance of each species
varN=function(v,nt,h,Xbar,N){
   out=lapply(1:Nspp,function(i,v,nt,h,Xbar,N) h[i]*sum((exp(v[[i]]-Xbar[i])^2)*nt[[i]])/N[i],v=v,nt=nt,h=h,Xbar=Xbar,N=N)
   return(unlist(out))
}  
              
# Function to do an image plot of a matrix in the usual orientation, A(1,1) at top left  
matrix.image=function(x,y,A,col=topo.colors(100),...) {
	nx=length(x); ny=length(y); 
	x1=c(1.5*x[1]-0.5*x[2],1.5*x[nx]-0.5*x[nx-1]); 
	y1=c(1.5*y[1]-0.5*y[2],1.5*y[ny]-0.5*y[ny-1]); 
	image(list(x=x,y=y,z=t(A)),xlim=x1,ylim=rev(y1),col=col,bty="u",...);  
}

#============================================================================================#
# (IV) Calculate the equilibrium areas.
#============================================================================================# 

## initial population density vector
nt=v
for(i in 1:Nspp) nt[[i]][]=0.1
new.nt=nt

# set up matrix to record cover
covSave = matrix(NA,tlimit,Nspp)
covSave[1,]=sumCover(v,nt,h,A)

# set up list to store size distributions
sizeSave=list(NULL)
for(i in 1:Nspp){
  sizeSave[[i]]=matrix(NA,length(v[[i]]),(tlimit))
  sizeSave[[i]][,1]=nt[[i]]/sum(nt[[i]])
}

# initial densities 
Nsave=matrix(NA,tlimit,Nspp)
Nsave[1,]=sumN(nt,h)

if(exists("climYrSave")){
 climYrSave=as.numeric(climYrSave[,1])
 if(length(climYrSave)!=tlimit) stop("Check length climYrSave")
}else{
 climYrSave=sample(1:Nyrs,tlimit,replace=T)
}

if(exists("randYrSave")){
 randYrSave=as.numeric(randYrSave[,1])
 if(length(randYrSave)!=tlimit) stop("Check length randYrSave")
}else{
 randYrSave=sample(1:Nyrs,tlimit,replace=T)
}

for (i in 2:(tlimit)){

  climYear=climYrSave[i]
  randYear=randYrSave[i]      # decouple climate effects and random year effects
 
  #get recruits per area
  cover=covSave[i-1,] ; N=Nsave[i-1,]
  recWeather=(climD[climYear,]-climDmean)/climDsd    # standardize
  rpa=get.rpa(Rpars,cover,weather=recWeather,doYear=randYear)
  
  #calculate size-specific crowding
  alpha=Spars$alpha 
  for(ii in 1:Nspp){ 
 	  # now do conspecific no overlap W
 	 	Ctot[ii]=h[ii]*sum(expv[[ii]]*nt[[ii]]) 
    Cr[[ii]]=splinefun(b.r[[ii]],h[ii]*c(0,cumsum(expv[[ii]]*nt[[ii]])),method="natural")
    Wmat[[ii]][,ii]=Wrij(v.r[[ii]],ii,ii)/A
  }
 
  for(doSpp in 1:Nspp){  
    if(cover[doSpp]>0){    
      # make kernels and project
      K.matrix=make.K.matrix(v[[doSpp]],Wmat[[doSpp]],Rpars,rpa,Gpars,Spars,weather=climD[climYear,],doYear=randYear,doSpp)	
      new.nt[[doSpp]]=K.matrix%*%nt[[doSpp]] 
      sizeSave[[doSpp]][,i]=new.nt[[doSpp]]/sum(new.nt[[doSpp]])  
    }    
  } # next species
	
  nt=new.nt 
	covSave[i,]=sumCover(v,nt,h,A)  # store the cover as cm^2/cm^2
  Nsave[i,]=sumN(nt,h)
 	print(i);flush.console()
 	if(sum(is.na(nt))>0) browser()  
 	
} # next time step

## Figures ==============================================================

par(mfrow=c(2,2),tcl=-0.2,mgp=c(2,0.5,0)) 
myCol=c("black","gold","blue","red")
#cover
boxplot(as.data.frame(100*covSave[(burn.in+1):tlimit,]),ylab="Cover (%)",names=sppList,col=myCol)
abline(h=0)
#density
boxplot(as.data.frame(Nsave[(burn.in+1):tlimit,]),ylab="Density",names=sppList,col=myCol)
abline(h=0)
# average size distribution  
plot(1,1,type="n",xlim=c(log(0.15),log(max(maxSize))),ylim=c(0,0.1),
  xlab="Size",ylab="Frequency")
for(i in 1:Nspp){
 lines(v[[i]],rowMeans(sizeSave[[i]][,(burn.in+1):tlimit]),col=myCol[i])
}
# example time series
matplot((burn.in+1):tlimit,100*covSave[(burn.in+1):tlimit,],type="l",col=myCol,
  xlab="Time",ylab="Cover (%)")

# plot changes in size distribution
X11()
par(mfrow=c(2,2),tcl=-0.2)
for(i in 1:Nspp) {image(sizeSave[[i]][],x=v[[i]],y=1:tlimit,
  xlab="size",ylab="time",main=sppList[i]) }

## Write ouput==============================================================

# output=data.frame(cbind(c(1:NROW(covSave)),covSave))
# names(output)=c("time",sppList)
# write.table(output,outfile1,row.names=F,sep=",")
# 
# for(i in 1:Nspp){
#  filename=paste(sppList[i],"_",outfile2,sep="")
#  tmp=rowMeans(sizeSave[[i]][,(burn.in+1):tlimit])
#  output=data.frame(v[[i]],tmp)
#  names(output)=c("size","freq")
#  write.table(output,filename,row.names=F,sep=",")
# }
# 
# #write years used
# write.table(climYrSave,"climYears.csv",row.names=F,sep=",")
# write.table(randYrSave,"randYears.csv",row.names=F,sep=",")
