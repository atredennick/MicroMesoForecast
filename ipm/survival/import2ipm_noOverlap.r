# import & format survival parameters
# then define survival function

climNames=names(climD)
climNames2=paste("logarea.",climNames,sep="")

# survival parameters
Spars=list(intcpt=rep(NA,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),intcpt.gr=matrix(0,6,Nspp),
  slope=rep(NA,Nspp),slope.yr=matrix(0,Nyrs,Nspp),
  clim=matrix(0,Nspp,length(climNames)),slopeXclim=matrix(0,Nspp,length(climNames)),
  nb=rep(0,Nspp),slopeXnb=rep(0,Nspp),alpha=matrix(NA,Nspp,Nspp))
colnames(Spars$clim)=climNames; colnames(Spars$slopeXclim)=climNames2

for(i in 1:Nspp){
  infile=paste("survival/Surv_params_",sppList[i],".csv",sep="")
  Sdata=read.csv(infile)
  Spars$intcpt[i]=Sdata$X.Intercept.[1]    # grand intercept
  Spars$intcpt.yr[,i]=Sdata$X.Intercept..yr # random effects on intercept
  tmp=which(names(Sdata)=="Group")    # group effects
  if(length(tmp)>0) Spars$intcpt.gr[,i]=Sdata$Group[!is.na(Sdata$Group)] # get spatial average
  Spars$slope[i]=Sdata$logarea[1]  #slope
  # random effects on slope
  tmp=which(names(Sdata)=="logarea.yr")
  if(length(tmp)>0) Spars$slope.yr[,i]=Sdata[,tmp]
  # climate effects
  tmp1=match(climNames,names(Sdata)) ; tmp1=tmp1[!is.na(tmp1)]
  tmp2=Sdata[1,tmp1]
  tmp3=match(names(Sdata)[tmp1],climNames)
  Spars$clim[i,c(tmp3)]=as.numeric(tmp2)
  # climate X slope interactions
  tmp1=match(climNames2,names(Sdata)) ; tmp1=tmp1[!is.na(tmp1)]
  tmp2=Sdata[1,tmp1]
  tmp3=match(names(Sdata)[tmp1],climNames2)
  Spars$slopeXclim[i,c(tmp3)]=as.numeric(tmp2)
  # get competition coefficients
  tmp=which(is.element(names(Sdata),"crowd"))
  if(length(tmp)>0) Spars$nb[i]=as.numeric(Sdata[1,tmp])
  # get competition X size interactions coefficients
  tmp=which(is.element(names(Sdata),"logarea.crowd"))
  if(length(tmp)>0) Spars$slopeXnb[i]=as.numeric(Sdata[1,tmp])   
  Spars$alpha[i,]=Sdata$alpha[1:length(sppList)]  #alphas
} # next i
rm(Sdata)

## survival function: probability an individual of size u survives  (u is on log scale)
S=function(u,W,Spars,weather,doYear,doSpp){
  mu=Spars$intcpt[doSpp]+Spars$intcpt.yr[doYear,doSpp]+(Spars$slope[doSpp]+Spars$slope.yr[doYear,doSpp])*u+
    sum(Spars$clim[doSpp,]*weather)+sum(Spars$slopeXclim[doSpp,]*weather)*u+
    W[,doSpp]*Spars$nb[doSpp]+(W[,doSpp]*u)*Spars$slopeXnb[doSpp]
  return(inv.logit(mu))
}


