#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

sppList=sort(c("BOGR","HECO","PASM","POSE"))
alpha.effect=c(0.010,0.025,0.020,0.036)  ##take this from multiple-species IPM
doSpp="HECO"
climD=read.csv("../../weather/Climate.csv")

outfile=paste("Growth_params_",doSpp,".csv",sep="")

growDfile=paste("../../speciesData/",doSpp,"/growDnoNA.csv",sep="")
growD=read.csv(growDfile)
growD$Group=as.factor(substr(growD$quad,1,1)) ##add Group information

D=growD  #subset(growD,allEdge==0)
D$logarea.t0=log(D$area.t0)
D$logarea.t1=log(D$area.t1)
D$quad=as.character(D$quad)
climD$year=climD$year-1900
D=merge(D,climD)
D$year=as.factor(D$year)


##then we moved some specific points:
tmp2<-which(D$quad=="A12" & D$year==44)
tmp3<-which(D$quad=="B1"  & D$year==44)
tmp41<-which(D$quad=="E4" & D$year==33) 
tmp42<-which(D$quad=="E4" & D$year==34) 
tmp43<-which(D$quad=="E4" & D$year==43)
tmp44<-which(D$quad=="E4" & D$year==44)

tmpONE<-c(tmp2,tmp3,tmp41,tmp42,tmp43,tmp44)
if(length(tmpONE)>0) D<-D[-tmpONE,]



# calculate crowding 
for(i in 1:length(sppList)){
  distDfile=paste("../../speciesData/",sppList[i],"/",sppList[i],"_genet_xy.csv",sep="")
  if(i==1){
    distD=read.csv(distDfile)
    distD$nbSpp=sppList[i]  
  }else{
    tmp=read.csv(distDfile)
    tmp$nbSpp=sppList[i] 
    distD=rbind(distD,tmp)
  }
}

distD=distD[,c("quad","year","trackID","area","nbSpp","x","y")]
W=matrix(NA,dim(D)[1],length(sppList))
for(i in 1:dim(D)[1]){
  tmpD=subset(distD,year==D$year[i] & quad==D$quad[i])
  focal=which(tmpD$trackID==D$trackID[i] & tmpD$nbSpp==doSpp)
  xx=tmpD$x[focal] ; yy=tmpD$y[focal]
  tmpD$distance=sqrt((xx-tmpD$x)^2+(yy-tmpD$y)^2)
  tmpD=subset(tmpD,distance>0)
  if(dim(tmpD)[1]>0){
    for(k in 1:length(sppList)){
      sppI=which(tmpD$nbSpp==sppList[k])
      if(length(sppI)>0){
        W[i,k]=sum(exp(-1*alpha.effect[k]*tmpD$distance[sppI]^2)*tmpD$area[sppI])         
      }else{
        W[i,k]=0
      }
    }
  }else{
    W[i,]=0
  }   
}

crowd=W[,which(sppList==doSpp)]

# compare fixed effects models
out1=glm(logarea.t1~Group+logarea.t0*crowd+  
           pptLag+logarea.t0:pptLag+
           ppt1+logarea.t0:ppt1+
           TmeanSpr1+logarea.t0:TmeanSpr1+
           ppt2+logarea.t0:ppt2+
           TmeanSpr2+logarea.t0:TmeanSpr2+
           ppt1:TmeanSpr1+logarea.t0:ppt1:TmeanSpr1+
           ppt2:TmeanSpr2+logarea.t0:ppt2:TmeanSpr2,
         data=D)
print(paste("m1: ",AIC(out1),sep=""))


library(MASS)
out2=stepAIC(out1,scope=list(upper=~.,lower=~logarea.t0+crowd),trace=T)

##Based on stepwise AIC, it is a bit tricky choosing the final model...
# fit final mixed effect model
library(lme4)
out=lmer(logarea.t1~logarea.t0*crowd+  
           pptLag+logarea.t0:pptLag+
           ppt1+logarea.t0:ppt1+
           TmeanSpr1+logarea.t0:TmeanSpr1+
           ppt2+logarea.t0:ppt2+
           TmeanSpr2+logarea.t0:TmeanSpr2+
           ppt1:TmeanSpr1+logarea.t0:ppt1:TmeanSpr1+
           ppt2:TmeanSpr2+logarea.t0:ppt2:TmeanSpr2+
           (1|Group)+(logarea.t0|year),data=D)

#check correlations between random year effects and climate
cor(cbind(ranef(out)$year,climD))

library(INLA)
#Set up ID variables for INLA random effects
D$yearID <- D$year #for random year offset on intercept
# D$groupID <- as.numeric(D$Group)
formula <- logarea.t1 ~ logarea.t0*crowd+
  pptLag+logarea.t0:pptLag+
  ppt1+logarea.t0:ppt1+
  TmeanSpr1+logarea.t0:TmeanSpr1+
  ppt2+logarea.t0:ppt2+
  TmeanSpr2+logarea.t0:TmeanSpr2+
  ppt1:TmeanSpr1+logarea.t0:ppt1:TmeanSpr1+
  ppt2:TmeanSpr2+logarea.t0:ppt2:TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea.t0, model="iid", prior="normal",param=c(1,0.001))

outINLA <- inla(formula, data=D,
                family=c("gaussian"), verbose=TRUE,
                control.compute=list(dic=T,mlik=T),
                control.predictor = list(link = 1))
summary(outINLA)

#fit variance
x = outINLA$summary.fitted.values$mean #fitted values from INLA
y = (D$logarea.t1-outINLA$summary.fitted.values$mean)^2 #calculates the variance (residuals^2)
plot(x,y)
outVar=nls(y~a*exp(b*x),start=list(a=1,b=0)) #model the variance as function of fitted size

#Collect parameters
#random year and group effects
params <- as.data.frame(outINLA$summary.random$yearID[,1:2])
params <- cbind(params, outINLA$summary.random$year[,2])
names(params) <- c("Year", "Intercept.yr", "logarea.t0.yr")
tmp <- as.data.frame(outINLA$summary.random$Group[,2])
names(tmp) <- "Group"
tmp$GrpName <- outINLA$summary.random$Group[,1]
tmp[(NROW(tmp)+1):NROW(params),]=NA
params=cbind(params,tmp)
#fixed effects
fixed <- as.data.frame(outINLA$summary.fixed)[,1:2]
tmp=matrix(NA,dim(params)[1],nrow(fixed))
colnames(tmp)=rownames(fixed)
tmp[1,]=fixed[,1]
params=cbind(params,tmp)
params$alpha=NA; params$alpha[1:length(sppList)]=alpha.effect
#variance 
params$sigma.a=NA; params$sigma.a[1]=coef(outVar)[1] 
params$sigma.b=NA; params$sigma.b[1]=coef(outVar)[2]

# write output
write.table(params,outfile,row.names=F,sep=",")
