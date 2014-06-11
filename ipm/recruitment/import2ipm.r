# all climate variables used
climNames=names(climD)

# import parameters
# recruitment parameters
Rpars=list(intcpt.mu=rep(0,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),betaC=matrix(0,Nspp,length(climNames)),
  intcpt.gr=matrix(0,6,Nspp),g.tau=rep(NA,Nspp),
  dd=rep(NA,Nspp),theta=rep(NA,Nspp),sizeMean=rep(NA,Nspp),sizeVar=rep(NA,Nspp),
  recSizes=list(1))

 infile="recruitment/recruit_params_Ryr.csv"    
 Rdata=read.csv(infile)
 # subset out non-essential parameters
 tmp=c(grep("lambda",row.names(Rdata)),grep("deviance",row.names(Rdata)),grep("DIC",row.names(Rdata)))   #group stuff?
 Rdata=Rdata[-tmp,]
 tmp=paste("Rpars$",row.names(Rdata),"<-",Rdata[,1],sep="")
 eval(parse(n=dim(Rdata)[1],text=tmp))
 for(i in 1:Nspp){
   infile=paste("H:/idahochart/ipm/speciesData/",sppList[i],"/recSize.csv",sep="")
   recSize=read.csv(infile)
   Rpars$sizeMean[i]=mean(log(recSize$area))
   Rpars$sizeVar[i]=var(log(recSize$area))
   #Rpars$recSizes[[i]]=recSize$area
 }
Rpars$dd=t(Rpars$dd) # c[i,j] = effect of j on i
rm(Rdata)

# define recruitment function
#number of recruits per area produced 
# cover is stored in absolute area (cm^2)
get.rpa=function(Rpars,cover,weather,doYear){
    # cover is in m^2 per m^2; convert to % scale:
    cover2=cover*100
    # calculate recruits
    Nspp=length(cover)
    mu=rep(NA,Nspp)
    for(i in 1:Nspp){
      mu[i]=cover2[i]*exp(Rpars$intcpt.yr[doYear,i]+sum(Rpars$betaC[i,]*weather)+sqrt(cover2[i])*Rpars$dd[i]) 
    }
    if(sum(is.na(mu))>0) browser() # stop for errors
    rpa=mu/(cover*A)  # convert from number recruits to recruits per cm^2
    return(rpa)
}

# Fecundity function, expected number of recruits of size y produced by a size x individual
# The size distribution of recruits is on the log scale
f=function(v,u,Rpars,rpa,doSpp) { 
  	nRecruits = rpa[doSpp]*exp(u)	
	  #probability of producing a seedling of size v
	  tmp=dnorm(v,Rpars$sizeMean[doSpp],sqrt(Rpars$sizeVar[doSpp]))/(1-pnorm(-1.61,Rpars$sizeMean[doSpp],sqrt(Rpars$sizeVar[doSpp])))
    #number recruits of each size 
    f=nRecruits*tmp
	  return(f)
}   

