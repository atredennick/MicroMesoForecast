# Individually-based model for a multiple species,
# with density-dependence and explicit space (random
# spatial pattern)

# Gaussian distance function

# This version uses the negative binomial recruitment function

# Initialized from observed size distributions

# PBA  5-21-2010

setwd("H:/idahochart/ipm/dd_climate_v4")
outfile1="equilibCover.csv"    
outfile2="sim_summary.csv"
outfile3="sim_xy.csv"
climateFile="Climate.csv"
totSims=2000
totT=100    # time steps of simulation
burn.in=20  # time steps to discard before calculating cover values
L=100   # dimension of square quadrat (cm)
expand=3    # 1 = 1x1 m^2, 2 = 2x2m^2, etc
init.cover=c(3,1,1,1)   # in % cover
maxSize=c(8000,500,500,500)
minSize=0.25
sppList=c("ARTR","HECO","POSE","PSSP")
Nyrs=22
myCol=c("black","gold","blue","red")
doGroup=NA  # NA for spatial avg., values 1-6 for a specific group

# GET VITAL RATE PARAMETERS & FUNCTIONS ------------------------------------------------
Nspp=length(sppList)

# get climate data
obsD=read.csv(obsClimateFile)
climD=obsD[,2:NCOL(obsD)]  #drop year column
climD$ppt1.TmeanSpr1=climD$ppt1*climD$TmeanSpr1
climD$ppt2.TmeanSpr2=climD$ppt2*climD$TmeanSpr2
climDmean=colMeans(climD); climDsd=apply(climD,MARGIN=2,FUN=sd)

# set up survival parameters and function
source("survival/import2ibm.r")
# set up growth parameters and function
source("growth/import2ibm.r")
# set up recruitment parameters and function
source("recruitment/import2ibm.r")

# model spatial group variation (or not)
if(!is.na(doGroup)){
  Spars$intcpt=Spars$intcpt+Spars$intcpt.gr[doGroup,]
  Gpars$intcpt=Gpars$intcpt+Gpars$intcpt.gr[doGroup,]
  Rpars$intcpt.yr=Rpars$intcpt.yr+matrix(Rpars$intcpt.gr[doGroup,],Nyrs,Nspp,byrow=T)
}

# get observed size distribution for initialization
size.init=list()
for(i in 1:Nspp){
  infile=paste("H:/idahochart/ipm/speciesData/",sppList[i],"/",sppList[i],"_genet_xy.csv",sep="")
  tmp=read.csv(infile)
  size.init[[i]]=tmp$area  
}

# OTHER FUNCTIONS---------------------------------------------------------
library(boot)
library(mvtnorm)
library(msm)  

# crowding function, assumes toroidal landscape
getDist=function(plants,L,expand){
 xdiff=abs(outer(plants[,3],plants[,3],FUN="-"))
 tmp=which(xdiff>((L*expand)/2))
 xdiff[tmp]=(L*expand)-xdiff[tmp]
 ydiff=abs(outer(plants[,4],plants[,4],FUN="-"))
 tmp=which(ydiff>((L*expand)/2))
 ydiff[tmp]=(L*expand)-ydiff[tmp]
 distMat=sqrt(xdiff^2+ydiff^2) 
 distMat[distMat==0]=NA
 return(distMat)
}

getCrowding=function(plants,alpha,distMat){
 if(dim(plants)[1]>1){  
   distMat=exp(-1*alpha[plants[,1]]*distMat^2)
   sizeMat=matrix(plants[,2],dim(plants)[1],dim(plants)[1])
   distSize=distMat*sizeMat
   out=sapply(1:4,function(i,distSize){ 
      colSums(matrix(distSize[plants[,1]==i,],sum(plants[,1]==i),NCOL(distSize)),na.rm=T)},
      distSize=distSize)
   out=t(out)
 }else{
   out=rep(0,Nspp)
 }
 out
}

# MAIN LOOP -------------------------------------------------------
outxy=matrix(NA,0,7)
colnames(outxy)=c("run","t","spp","size","x","y","id")
output=matrix(NA,0,3+2*Nspp)
colnames(output)=c("run","time","yrParams",paste("Cov.",sppList,sep=""),
  paste("N.",sppList,sep=""))  
for(iSim in 1:totSims){ 
  # INITIALIZE by drawing from observed sizes until target cover reached
  spp=NULL ; size=NULL
  for(iSpp in 1:Nspp){
    if(init.cover[iSpp]>0){
      n.init=round(init.cover[iSpp]*L/mean(size.init[[iSpp]])*expand^2)    
      target=init.cover[iSpp]*L*expand^2
      lower=target-0.1*target
      upper=target+0.1*target
      success=F
      while(success==F){
        sizeTry=sample(size.init[[iSpp]],n.init)
        if(sum(sizeTry)>lower & sum(sizeTry)<upper) success=T
      }
      size=c(size,sizeTry)
      spp=c(spp,rep(iSpp,length(sizeTry)))
    }
  }
  x=runif(length(size),0,L*expand) ; y=runif(length(size),0,L*expand)
  id=rep(1:length(size))
  plants=cbind(spp,size,x,y,id)
  lastID=max(plants[,5])
  
  # vectors for cover and density
  N=matrix(0,totT,Nspp)
  N=colSums(matrix(plants[,1],dim(plants)[1],Nspp)==matrix(1:Nspp,dim(plants)[1],Nspp,byrow=T))
  A=rep(0,Nspp)
  for(i in 1:Nspp){
    A[i]=sum(plants[plants[,1]==i,2])/(expand^2*L^2)
  }
  new.N=rep(0,Nspp); new.A=rep(0,Nspp)
  tmp=c(iSim,1,0,A,N)
  output=rbind(output,tmp)
  
  # plot initial conditions
  par(mgp=c(2,0.5,0),tcl=-0.2)
  symbols(x = plants[,3], y = plants[,4], circles = sqrt(plants[,2]/pi),fg=myCol[plants[,1]],
    xlim=c(0,L*expand),ylim=c(0,L*expand),main ="Time=1",xlab="x",ylab="y",inches=F,lwd=2)
  
  for(tt in 2:(totT)){
     
     # draw year effects
     climYr=sample(1:NROW(climD),1)
     doYr=sample(1:NROW(climD),1)
     nextplants=plants
     
     # recruitment
     recWeather=(climD[climYr,]-climDmean)/climDsd    # standardize
     newplants=recruit(Rpars,sizes=plants[,2],spp=plants[,1],weather=recWeather,doYr,lastID=lastID,L,expand)
            
     for(ss in 1:Nspp){
      if(N[ss]>0){ # make sure spp ss is not extinct
        
        # crowding
        distMat=getDist(plants,L,expand)
        W=getCrowding(plants,Gpars$alpha[ss],distMat)
        
        # growth
        newsizes=grow(Gpars,doSpp=ss,weather=climD[doYr,],doYr,sizes=plants[,2],crowding=W[ss,])
        if(is.na(sum(newsizes))) browser()
        if(sum(newsizes==Inf)>0) browser()
        
        # survival
        # uses same W as growth
        live=survive(Spars,doSpp=ss,weather=climD[doYr,],doYr,sizes=plants[,2],crowding=W[ss,])
  
        # combine growth and survival
        tmp=which(plants[,1]==ss)  # only alter plants of focal spp        
        nextplants[tmp,2]=newsizes[tmp]*live[tmp]   #update with G and S
        
       } # end if no plants
     } # next ss 
 
     nextplants=nextplants[nextplants[,2]>0,]    # remove dead plants 
     nextplants=rbind(nextplants,newplants)     # add recruits
     
     # output cover and density
     A[]=0; N[]=0
     tmp=aggregate(nextplants[,2],by=list(nextplants[,1]),FUN=sum)
     A[tmp[,1]]=tmp[,2]/(expand^2*L^2)
     tmp=aggregate(rep(1,dim(nextplants)[1]),by=list(nextplants[,1]),FUN=sum)
     N[tmp[,1]]=tmp[,2]/(expand^2)
    
     lastID=max(nextplants[,5])
     plants=nextplants
     if(is.matrix(plants)==F) plants=matrix(plants,nrow=1,ncol=length(plants))
     # plot
     if(sum(N)>0){
      symbols(x = plants[,3], y = plants[,4], circles = sqrt(plants[,2]/pi),fg=myCol[plants[,1]],
        xlim=c(0,L*expand),ylim=c(0,L*expand),main =paste("Time=",tt,sep=""),
        xlab="x",ylab="y",inches=F,lwd=2)
     }else{
       break
     }
     output=rbind(output,c(iSim,tt,doYr,A,N)) 
     
     # save xy coordinates for spatial analysis
    if(tt>burn.in & tt%%10==0){
      if(sum(plants[,1]==1)>4) {
        tmp=cbind(rep(iSim,dim(plants)[1]),rep(tt,dim(plants)[1]),plants)
        outxy=rbind(outxy,tmp)
      }
    }
      
  } # next tt
  print(paste("Sim ",iSim," complete",sep=""))
  print(date())
  flush.console()
} # next iSim

# write output
row.names(output)=1:dim(output)[1]
write.table(output,outfile1,row.names=F,sep=",")

row.names(outxy)=1:dim(outxy)[1]
write.table(outxy,outfile3,row.names=F,sep=",")

# get extinction times for each species in each run
extinct=as.data.frame(output[,1:(3+Nspp)])
extinct=subset(extinct,time>1) # get rid of initialization years
eTimes=data.frame("run"=1:totSims)
for(i in 1:Nspp){
 tmpD=extinct[,c(1,2,(3+i))]
 tmpD=subset(tmpD,tmpD[,3]==0)
 if(dim(tmpD)[1]>0){
   tmpD=aggregate(tmpD$time,by=list(run=tmpD$run),FUN=min)
   names(tmpD)[2]=sppList[i]
 }else{
   tmpD=tmpD[,1:2]
   tmpD[1,]=c(1,NA)
   names(tmpD)[2]=sppList[i]
 }
 eTimes=merge(eTimes,tmpD,all.x=T)
}
freq.extinct=colSums(!is.na(eTimes[,2:NCOL(eTimes)]),na.rm=T)/totSims
mean.eTimes=colMeans(eTimes[,2:NCOL(eTimes)],na.rm=T)
sd.eTimes=apply(eTimes[,2:NCOL(eTimes)],MARGIN=2,FUN=sd,na.rm=T)
rm(extinct)

# get average cover for non-extinction runs
cover=as.data.frame(output)
cover=subset(cover,time>burn.in)  # remove burn in period
#cover[cover>1]=NA  # remove outliers
index=which(rowSums(is.na(eTimes[,2:NCOL(eTimes)]))==Nspp)
index=which(is.element(cover$run,index))
cover=cover[index,]
mean.cover=colMeans(cover[,4:(3+Nspp)])
sd.cover=apply(cover[,4:(3+Nspp)],MARGIN=2,FUN=sd)

summary=data.frame(cbind(mean.cover,sd.cover,mean.eTimes,sd.eTimes,freq.extinct))
summary$totT=NA
summary$totT[1]=totT
summary$totRuns=NA
summary$totRuns[1]=totSims
write.table(summary,outfile2,row.names=F,sep=",")

# figures
par(mfrow=c(1,3),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,3,1,1))
barplot(mean.cover*100,names.arg=sppList,ylab="Mean cover (%)",
  ylim=c(0,max(mean.cover+sd.cover)*100))
xx=as.numeric(barplot(mean.cover*100,plot=F))
arrows(x0=xx,y0=(mean.cover*100),x1=xx,y1=((mean.cover+sd.cover)*100),
  code=2,angle=90,length=0.05)

barplot(freq.extinct,names.arg=sppList,ylim=c(0,1),
  ylab="Extinction frequency")

barplot(mean.eTimes,names.arg=sppList,ylab="Mean extinction time",
  ylim=c(0,totT))






