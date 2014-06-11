# Import and format growth parameters
# then define growth function
climNames=names(climD)
climNames2=paste("logarea.t0.",climNames,sep="")

# growth parameters
Gpars=list(intcpt=rep(NA,Nspp),intcpt.gr=matrix(0,6,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),
  slope=rep(NA,Nspp),slope.yr=matrix(0,Nyrs,Nspp),
  clim=matrix(0,Nspp,length(climNames)),slopeXclim=matrix(0,Nspp,length(climNames)),
  nb=rep(0,Nspp),alpha=matrix(NA,Nspp,Nspp),
  sigma2.a=rep(NA,Nspp),sigma2.b=rep(NA,Nspp))
colnames(Gpars$clim)=climNames

for(i in 1:Nspp){
  infile=paste("growth/Growth_params_",sppList[i],".csv",sep="")
  Gdata=read.csv(infile)
  Gpars$intcpt[i]=Gdata$X.Intercept.[1]    # grand intercept
  Gpars$intcpt.yr[,i]=Gdata$X.Intercept..yr       # random year effects
  tmp=which(names(Gdata)=="Group")    # group effects
  if(length(tmp)>0) Gpars$intcpt.gr[,i]=Gdata$Group[!is.na(Gdata$Group)] # get spatial average
  Gpars$slope[i]=Gdata$logarea.t0[1]  #slope
  # random effects on slope
  tmp=which(names(Gdata)=="logarea.t0.yr")
  if(length(tmp)>0) Gpars$slope.yr[,i]=Gdata[,tmp]
  # climate effects
  tmp1=match(climNames,names(Gdata)) ; tmp1=tmp1[!is.na(tmp1)]
  tmp2=Gdata[1,tmp1]
  tmp3=match(names(Gdata)[tmp1],climNames)
  Gpars$clim[i,c(tmp3)]=as.numeric(tmp2)
  # climate X slope interactions
  tmp1=match(climNames2,names(Gdata)) ; tmp1=tmp1[!is.na(tmp1)]
  tmp2=Gdata[1,tmp1]
  tmp3=match(names(Gdata)[tmp1],climNames2)
  Gpars$slopeXclim[i,c(tmp3)]=as.numeric(tmp2)
  # get competition coefficients
  tmp=which(is.element(names(Gdata),"crowd"))
  if(length(tmp)>0) Gpars$nb[i]=as.numeric(Gdata[1,tmp])
  Gpars$alpha[i,]=Gdata$alpha[1:length(sppList)]  #alphas
  Gpars$sigma2.a[i]=Gdata$sigma.a[1]
  Gpars$sigma2.b[i]=Gdata$sigma.b[1]
} # next i
rm(Gdata)

# growth function
G=function(v,u,W,Gpars,weather,doYear,doSpp){
  mu=Gpars$intcpt[doSpp]+Gpars$intcpt.yr[doYear,doSpp]+(Gpars$slope[doSpp]+Gpars$slope.yr[doYear,doSpp])*u+
    sum(Gpars$clim[doSpp,]*weather)+sum(Gpars$slopeXclim[doSpp,]*weather)*u+
    Gpars$nb[doSpp]*W[,doSpp]
  sigma2=Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  out=dnorm(v,mu,sqrt(sigma2))
  out
}


