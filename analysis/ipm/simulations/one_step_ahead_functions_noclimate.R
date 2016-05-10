
#============================================================================================#
# Utility functions
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
                       Rpars,rpa,Gpars,Spars,doYear,doSpp)  #vital rate arguments
{
  f(v,u,Rpars,rpa,doSpp)+S(u,muWS,Spars,doYear,doSpp)*G(v,u,muWG,Gpars,doYear,doSpp)
}

# Function to make iteration matrix based only on mean crowding
make.K.matrix=function(v,muWG,muWS,Rpars,rpa,Gpars,Spars,doYear,doSpp) {
  #muWS=expandW(v,v,muWS)  # for multispecies models
  #muWG=expandW(v,v,muWG)
  muWS<-rep(muWS,length(v))
  muWG<-rep(muWG,length(v))
  K.matrix=outer(v,v,make.K.values,muWG,muWS,Rpars,rpa,Gpars,Spars,doYear,doSpp)
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

# build kernel and project population
projectIPM<-function(nt,doYear,doGroup,sppCode){
  
  # get vital rate parameters
  Spars<-getSurvCoefs(doYear,doGroup)
  Gpars<-getGrowCoefs(doYear,doGroup)
  Rpars<-getRecCoefs(doYear,doGroup)
  Rpars$sizeMean <- rec_size_mean
  Rpars$sizeVar <- rec_size_var
  
  # calculate intraspecific crowding 
  #Ctot=h*sum(expv*nt) #total cover
 # Cr=splinefun(b.r,h*c(0,cumsum(expv*nt)),method="natural") #Cr is a function      
  WfunG=splinefun(size.range,WriG(size.range))
  WmatG=WfunG(v.r)/Atotal
  WfunS=splinefun(size.range,WriS(size.range))
  WmatS=WfunS(v.r)/Atotal
  
  # get recruits per area
  cover.lag<-sum(exp(v)*nt)/Atotal  # sumCover(v,nt,h,Atotal) multiply by h or not???
  rpa<-rep(0,n_spp); rpa[sppCode]=get_rpa(Rpars,cover.lag)  # set up this way to eventually run multiple species at once
  
  # make kernel and project population
  K.matrix=make.K.matrix(v,WmatG,WmatS,Rpars,rpa,Gpars,Spars,doYear,sppCode)  
  new.nt=K.matrix%*%nt
  
  return(new.nt)
  
} 

