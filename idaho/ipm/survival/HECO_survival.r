# PBA 6-9-14

setwd("H:\\idahochart\\ipm\\dd_climate_v4\\survival")
dataDir="H:\\idahochart\\ipm\\speciesData\\"
sppList=(c("ARTR","HECO","POSE","PSSP"))
doSpp="HECO"
alpha.effect=c(0.012,0.06,0.05,0.03) # for spp in alphabetical order
climD=read.csv("H:\\idahochart\\ipm\\dd_climate_v4\\Climate.csv")
#--------------------------------------

outfile=paste("Surv_params_",doSpp,".csv",sep="")

survDfile=paste(dataDir,doSpp,"\\survD.csv",sep="")
survD=read.csv(survDfile)
D=survD #subset(survD,allEdge==0)
D$logarea=log(D$area)
D$quad=as.character(D$quad)
climD$year=climD$year-1900
D=merge(D,climD)
D$year=as.factor(D$year)

# calculate crowding
ringD=read.csv(paste("H:\\idahochart\\ipm\\speciesData\\",doSpp,"\\",doSpp,"_nbhood_rings.csv",sep=""))
ringD$year<-as.factor(ringD$year)
midRings <- seq(2.5,142.5,5)  # could be extracted from ringD names, hardwired here for convenience
Nrings <- length(midRings)
D<-merge(D,ringD,by.x=c("quad","year","trackID"),by.y=c("quad","year","genetID"))
D=D[order(D$X),]
dist.wts=exp(-alpha.effect[which(sppList==doSpp)]*midRings^2) 
sppCols=which(substr(names(D),1,4)==doSpp); 
sppData=as.matrix(D[,sppCols]); 
crowd=sppData%*%dist.wts
 
#final mixed effect model
library(lme4)

# no Group effect, use climate variables from out2
out=glmer(survives~logarea+crowd+
    pptLag+ppt1+TmeanSpr1+ ppt2+TmeanSpr2+
    (crowd+logarea|year)+(1|Group),data=D,family=binomial,
    verbose=TRUE, control=list(msVerbose=TRUE, maxIter=1000, maxFN=2500)) 


# year random effects
params=as.data.frame(ranef(out)$year)
names(params)=paste(names(params),".yr",sep="")
#group effects
tmp=as.data.frame(ranef(out)$Group)
names(tmp)="Group"
tmp$GrpName=row.names(tmp)
tmp[(NROW(tmp)+1):NROW(params),]=NA
params=cbind(params,tmp)
#fixed effects
tmp=matrix(NA,dim(params)[1],length(fixef(out)))
colnames(tmp)=names(fixef(out))
tmp[1,]=fixef(out)
params=cbind(params,tmp)
params$alpha=NA; params$alpha[1:length(sppList)]=alpha.effect
write.table(params,outfile,row.names=F,sep=",")
