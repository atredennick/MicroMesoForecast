# try discrete distributions for number recruits produced

setwd("H:\\idahochart\\ipm\\dd_climate_v4\\recruitment\\")
sppList=c("ARTR","HECO","POSE","PSSP")
outfile="recruit_params_RYr.csv"
climD=read.csv("H:\\idahochart\\ipm\\dd_climate_v4\\Climate.csv")
climVars=c("pptLag","ppt1","ppt2","TmeanSpr1","TmeanSpr2")  
#--------------------------------------------------------

# subset climate data
tmp=which(is.element(names(climD),climVars))
climD=climD[,tmp]
#climD$pptXtmp1=climD$ppt1*climD$TmeanSpr1
#climD$pptXtmp2=climD$ppt2*climD$TmeanSpr2
# standardize climD
tmp1=matrix(colMeans(climD),NROW(climD),NCOL(climD),byrow=T)
tmp2=matrix(apply(climD,MARGIN=2,FUN=sd),NROW(climD),NCOL(climD),byrow=T)
climD=(climD-tmp1)/tmp2

# get recruitment data
Nspp=length(sppList)
for(i in 1:Nspp){
  infile1=paste("H:\\idahochart\\ipm\\speciesData\\",sppList[i],"\\recArea.csv",sep="")
  tmpD=read.csv(infile1)
  tmpD=tmpD[,c("quad","year","NRquad","totParea","Group")]
  names(tmpD)[3]=paste("R.",sppList[i],sep="")
  names(tmpD)[4]=paste("cov.",sppList[i],sep="")
  if(i==1){
    D=tmpD
  }else{
    D=merge(D,tmpD,all=T)
  }
}
D[is.na(D)]=0  # replace missing values 


# calculate mean cover by group and year
tmpD=D[,c("quad","year","Group",paste("cov.",sppList,sep=""))]
tmpD=aggregate(tmpD[,4:NCOL(tmpD)],by=list("year"=tmpD$year,
  "Group"=tmpD$Group),FUN=mean)
names(tmpD)[3:NCOL(tmpD)]=paste("Gcov.",sppList,sep="")

D=merge(D,tmpD,all.x=T)
  
y=as.matrix(D[,c(paste("R.",sppList,sep=""))])
R.tot=rowSums(y)
parents1=as.matrix(D[,c(paste("cov.",sppList,sep=""))])/100
parents2=as.matrix(D[,c(paste("Gcov.",sppList,sep=""))])/100
year=as.numeric(as.factor(D$year))
Nyears=length(unique(D$year))
N=dim(D)[1]
Nspp=length(sppList)
Group=as.numeric(D$Group)
Ngroups=length(unique(Group))
#weather=as.matrix(climD[,2:NCOL(climD)])
weather=as.matrix(climD)
tmp=matrix(colMeans(weather),NROW(weather),NCOL(weather),byrow=T)
weather=weather-tmp # deviations

if(dim(weather)[1]!=Nyears) stop("Check number of years & climate data")
Nclim=dim(weather)[2]

# plots
pdf("recruit_data.pdf",height=6,width=8)
par(mfrow=c(1,2),tcl=-0.2,mgp=c(2,0.5,0),mar=c(3,3,3,1))
wts=c(0.6,1,0.65,0.9)
for(i in 1:Nspp){
 plot(parents1[,i],y[,i],xlab="Local parents (% cover)",ylab="Recruits",main=sppList[i],pch=year,col=year)
 trueparents=wts[1]*parents1[,i]+(1-wts[1])*parents2[,i]
 plot(trueparents,y[,i],xlab="Mixed parents (% cover)",ylab="Recruits",main=sppList[i],pch=year,col=year)
}
dev.off()


# fit as negative binomial with random effects in WinBUGS
library(boot)
library(R2WinBUGS)

data=list("N","y","parents1","parents2",
  "year","Nyears","Nspp","Ngroups","Group","weather","Nclim")

inits=list(1)
inits[[1]]=list(betaC=matrix(-0.01,Nspp,Nclim),
  intcpt.mu=rep(1,Nspp),intcpt.yr=matrix(-1,Nyears,Nspp),yr.tau=rep(1,Nspp),
  intcpt.gr=matrix(-1,Ngroups,Nspp),g.tau=rep(1,Nspp),dd=rep(0,Nspp),theta=rep(1,Nspp)) 
inits[[2]]=list(betaC=matrix(0,Nspp,Nclim),
  intcpt.mu=rep(0,Nspp),intcpt.yr=matrix(0,Nyears,Nspp),yr.tau=rep(0.1,Nspp),
  betaC=matrix(0,Nspp,Nclim),g.tau=rep(0.1,Nspp),dd=rep(0,Nspp),theta=rep(2,Nspp))
  
params=c("betaC","intcpt.mu","intcpt.yr","yr.tau",
  "intcpt.gr","g.tau","dd","theta","u","lambda") 

out=bugs(data,inits,params,
  model.file="bugsRYr.txt",
  n.chains=2,
  n.iter=100000,
  n.burnin=20000,
  n.thin=100, 
  debug=T,DIC=T,
  bugs.directory="C:/WinBUGS14/")  
  
tmp=grep("lambda",row.names(out$summary))
A=row.names(out$summary)[tmp]
B=out$summary[tmp,1]
lambda=matrix(NA,dim(y)[1],Nspp)
C=paste(A,"<-",B,sep="")
eval(parse(n=length(A),text=C))
lambda[is.na(lambda)]=0
par(mfrow=c(2,2))
for(i in 1:Nspp){
  plot(y[,i],lambda[,i],xlab="Obs",ylab="Pred",main=sppList[i])
}
par(mfrow=c(2,2))
for(i in 1:Nspp){
  plot(parents1[,i],lambda[,i],xlab="Parents",ylab="Pred",main=sppList[i])
}

write.table(out$summary,outfile,row.names=T,sep=",")
tmp=paste("DIC",out$DIC,sep=",")
write.table(tmp,outfile,col.names=F,row.names=F,append=T)
