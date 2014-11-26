# try discrete distributions for number recruits produced
#clear everything, just to be safe 
rm(list=ls(all=TRUE))

sppList=sort(c("BOGR","HECO","PASM","POSE"))
outfile="recruit_params_GnoG.csv"


# get recruitment data
Nspp=length(sppList)
for(i in 1:Nspp){
  infile1=paste("../../../speciesData/",sppList[i],"/recArea.csv",sep="")
  tmpD=read.csv(infile1)
  
  tmpD$Group=as.factor(substr(tmpD$quad,1,1)) #add by Chengjin

  
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


##then we moved some specific points:
# tmp2<-which(D$quad=="A12" & D$year==44)
# tmp3<-which(D$quad=="B1"  & D$year==44)
# 
# 
# tmpONE<-c(tmp2,tmp3)
# if(length(tmpONE)>0) D<-D[-tmpONE,]


# calculate mean cover by group and year
tmpD=D[,c("quad","year","Group",paste("cov.",sppList,sep=""))]
tmpD=aggregate(tmpD[,4:NCOL(tmpD)],by=list("year"=tmpD$year,"Group"=tmpD$Group),FUN=mean)
names(tmpD)[3:NCOL(tmpD)]=paste("Gcov.",sppList,sep="")

D=merge(D,tmpD,all.x=T)




###if using square transform, there is no need of the following code
####################################################################
###here you need to check: both parents1 and parents2 are equal to 0 at the same time
parents1=as.matrix(D[,c(paste("cov.",sppList,sep=""))])/100 ##convert from absolute cover to [1,100] range
parents2=as.matrix(D[,c(paste("Gcov.",sppList,sep=""))])/100

##for species 1
tmp1L=which(parents1[,1]==0) ##lcoal
tmp1G=which(parents2[,1]==0) ##Group
tmp1=intersect(tmp1L,tmp1G)
##for species 2
tmp2L=which(parents1[,2]==0)
tmp2G=which(parents2[,2]==0)
tmp2=intersect(tmp2L,tmp2G)
##for species 3
tmp3L=which(parents1[,3]==0)
tmp3G=which(parents2[,3]==0)
tmp3=intersect(tmp3L,tmp3G)
##for species 4
tmp4L=which(parents1[,4]==0)
tmp4G=which(parents2[,4]==0)
tmp4=intersect(tmp4L,tmp4G)

tmp<-unique(c(tmp1,tmp2,tmp3,tmp4))
 
if(length(tmp)>0){
  parents1<-parents1[-tmp,] ##remove them
  parents2<-parents2[-tmp,] ##remove them
  y=as.matrix(D[,c(paste("R.",sppList,sep=""))])[-tmp,] ##remove them  
  year=as.numeric(as.factor(D$year))[-tmp] ##remove them
  Nyrs=length(unique(D$year))
  N=dim(D)[1]-length(tmp) ##reduce
  Nspp=length(sppList)
  Group=as.numeric(as.factor(D$Group))[-tmp] ##remove them ##first turn it as FACTOR, then to NUMERIC
  Ngroups=length(unique(Group))
} else {
  y=as.matrix(D[,c(paste("R.",sppList,sep=""))])
  year=as.numeric(as.factor(D$year))
  Nyrs=length(unique(D$year))
  N=dim(D)[1]
  Nspp=length(sppList)
  Group=as.numeric(as.factor(D$Group)) ##first turn it as FACTOR, then to NUMERIC
  Ngroups=length(unique(Group))
}


# fit as negative binomial with random effects in WinBUGS
library(coda)
library(rjags)

dataJ = list(N=N, 
             Y=y, 
             parents1=parents1, 
             parents2=parents2, 
             yrs=year, 
             nYrs=Nyrs,
             nSpp=Nspp, 
             nGrp=Ngroups, 
             grp=Group)

inits=NULL
inits[[1]]=list(intcpt.yr=matrix(0.1,Nyrs,Nspp),intcpt.mu=rep(0.5,Nspp),intcpt.tau=rep(1,Nspp),
  intcpt.gr=matrix(1,Ngroups,Nspp),g.tau=rep(1,Nspp),
  dd=matrix(-0.1,Nspp,Nspp),theta=rep(2,Nspp)) 
inits[[2]]=list(intcpt.yr=matrix(0.5,Nyrs,Nspp),intcpt.mu=rep(0.2,Nspp),intcpt.tau=rep(10,Nspp),
  intcpt.gr=matrix(0,Ngroups,Nspp),g.tau=rep(0.1,Nspp),
  dd=matrix(-0.05,Nspp,Nspp),theta=rep(2,Nspp))
  

params=c("intcpt.yr","intcpt.mu","intcpt.tau","intcpt.gr","g.tau","dd","theta","u","lambda") 

modelFile <- "recruitJAGS.R"

n.Adapt <- 1000
n.Up <- 1000
n.Samp <- 2000

jm <- jags.model(modelFile, data=dataJ, n.chains=length(inits),inits = inits, n.adapt = n.Adapt)
update(jm, n.iter=n.Up)
# Get sums of square residuals
zm <- coda.samples(jm, variable.names=params, n.iter=n.Samp, n.thin=10)
  
zmStat <- summary(zm)$stat
tmp=grep("lambda",row.names(zmStat))
A=row.names(zmStat)[tmp]
B=zmStat[tmp,1]
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

write.csv(zmStat,outfile,row.names=T)
# tmp=paste("DIC",out$DIC,sep=",")
# write.table(tmp,outfile,col.names=F,row.names=F,append=T)
