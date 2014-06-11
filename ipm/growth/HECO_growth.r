
setwd("H:\\idahochart\\ipm\\dd_climate_v4\\growth")
dataDir="H:\\idahochart\\ipm\\speciesData\\"
sppList=sort(c("PSSP","HECO","POSE","ARTR"))
alpha.effect=c(0.012,0.06,0.05,0.03) # for spp in alphabetical order
doSpp="HECO"
climD=read.csv("H:\\idahochart\\ipm\\dd_climate_v4\\Climate.csv")
#--------------------------------------

outfile=paste("Growth_params_",doSpp,".csv",sep="")

growDfile=paste(dataDir,doSpp,"\\growDnoNA.csv",sep="")
growD=read.csv(growDfile)
D=growD  #subset(growD,allEdge==0)
D$logarea.t0=log(D$area.t0)
D$logarea.t1=log(D$area.t1)
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


# fit final mixed effect model
library(lme4)

# model 7 is best, but behaves poorly, use model 5
out=lmer(logarea.t1~logarea.t0+crowd+
        pptLag+ppt1+TmeanSpr1+ppt2+TmeanSpr2+   
        (logarea.t0|year)+(1|Group),D)     

#check correlations between random year effects and climate
cor(cbind(ranef(out)$year,climD))

# fit variance
x=fitted(out)
y=resid(out)^2
plot(x,y)
outVar=nls(y~a*exp(b*x),start=list(a=1,b=0))

# write output
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
#variance 
params$sigma.a=NA; params$sigma.a[1]=coef(outVar)[1] 
params$sigma.b=NA; params$sigma.b[1]=coef(outVar)[2]
write.table(params,outfile,row.names=F,sep=",")
