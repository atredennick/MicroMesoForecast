# try discrete distributions for number recruits produced
#clear everything, just to be safe 
rm(list=ls(all=TRUE))

library(plyr)
library(reshape2)
library(coda)
library(rjags)
load.module("dic")

sppList=sort(c("BOGR","HECO","PASM","POSE"))

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
  year=D$year[-tmp] ##remove them
  Nyrs=length(unique(D$year))
  N=dim(D)[1]-length(tmp) ##reduce
  Nspp=length(sppList)
  Group=as.numeric(as.factor(D$Group))[-tmp] ##remove them ##first turn it as FACTOR, then to NUMERIC
  Ngroups=length(unique(Group))
} else {
  y=as.matrix(D[,c(paste("R.",sppList,sep=""))])
  year=D$year
  Nyrs=length(unique(D$year))
  N=dim(D)[1]
  Nspp=length(sppList)
  Group=as.numeric(as.factor(D$Group)) ##first turn it as FACTOR, then to NUMERIC
  Ngroups=length(unique(Group))
}

tmpY <- melt(y)
tmpP1 <- melt(parents1)
tmpP2 <- melt(parents2)
allD <- data.frame(species=tmpY$Var2,
                   year=rep(year,4),
                   group=rep(Group,4),
                   recruits=tmpY$value,
                   parents1=tmpP1$value,
                   parents2=tmpP2$value)

climD <- read.csv("../../../weather/Climate.csv")
climD[3:6] <- scale(climD[3:6], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900

allD <- merge(allD, climD)


# fit as negative binomial with random effects in JAGS
dataJ = list(nObs=nrow(allD), 
             Y=allD$recruits, 
             parents1=allD$parents1, 
             parents2=allD$parents2, 
             yrs=as.numeric(as.factor(allD$year)), 
             nYrs=Nyrs,
             nSpp=Nspp, 
             nGrp=Ngroups, 
             spp=as.numeric(as.factor(allD$species)),
             grp=allD$group,
             TmeanSpr1 = allD$TmeanSpr1,
             TmeanSpr2 = allD$TmeanSpr2,
             ppt1 = allD$ppt1,
             ppt2 = allD$ppt2)

iterations <- 50000
adapt <- 10000
mod <- jags.model("recruitmentAllSppCLIMATE_JAGS.R", data=dataJ, n.chains=3, n.adapt = adapt)
update(mod, n.iter=iterations)
dic <- jags.samples(mod, c("deviance"),
                    n.iter=iterations, n.thin=10)

write.csv(outDeviance, file="recruitmentDevianceCLIMATE.csv")



