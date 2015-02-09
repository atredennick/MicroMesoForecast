
# single species IPM with climate covariates
# this script simulates equilibrium cover

# simulates one time step for a give quad, year, and initial conditions


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







for (i in 2:(tlimit)){
  
  doYear<-yrSave[i]
  weather<-clim_data[clim_data$year==(1900+doYear),2:5]
  mcDraw<-sample(1:nMCMC,1)
  
  # get vital rate parameters
  Spars<-getSurvCoefs(doYear,mcDraw,doGroup)
  Gpars<-getGrowCoefs(doYear,mcDraw,doGroup)
  Rpars<-getRecCoefs(doYear,mcDraw,doGroup)
  Rpars$sizeMean <- rec_size_mean
  Rpars$sizeVar <- rec_size_var
  
  # calculate intraspecific crowding 
  Ctot=h*sum(expv*nt) #total cover
  Cr=splinefun(b.r,h*c(0,cumsum(expv*nt)),method="natural") #Cr is a function      
  WfunG=splinefun(size.range,WriG(size.range))
  WmatG=WfunG(v.r)/Atotal
  WfunS=splinefun(size.range,WriS(size.range))
  WmatS=WfunS(v.r)/Atotal
  
  # get recruits per area
  cover<-covSave[i-1]
  rpa<-rep(0,n_spp); rpa[sppCode]=get_rpa(Rpars,cover,weather)  # set up this way to eventually run multiple species at once
  
  # make kernel and project population
  K.matrix=make.K.matrix(v,WmatG,WmatS,Rpars,rpa,Gpars,Spars,doYear,sppCode,weather)  
  new.nt=K.matrix%*%nt
  sizeSave[,i]=new.nt/sum(new.nt)  
  
  # update and store state variables
  nt=new.nt
  covSave[i]<-sumCover(v,nt,h,Atotal)  # store the cover as cm^2/cm^2
  Nsave[i]=sumN(nt,h)
  
  if(sum(is.na(nt))>0) browser() # check for errors

  print(i)# print the simulation steps
  flush.console()
  
} # next time step



