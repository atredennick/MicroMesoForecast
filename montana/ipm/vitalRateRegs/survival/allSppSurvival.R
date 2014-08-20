#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

sppList=sort(c("BOGR","HECO","PASM","POSE"))
alpha.effect=c(0.010,0.025,0.020,0.036)  ##take this from multiple-species IPM

fixedEffects <- list()
SSRs <- matrix(ncol=5, nrow=(length(sppList)))

for(spp in 1:length(sppList)){
# for(spp in 2:2){
  doSpp=sppList[spp]
  climD=read.csv("../../../weather/Climate.csv")
  
  outfile=paste("Surv_params_",doSpp,".csv",sep="")
  
  growDfile=paste("../../../speciesData/",doSpp,"/survD.csv",sep="")
  growD=read.csv(growDfile)
  growD$Group=as.factor(substr(growD$quad,1,1)) ##add Group information
  
  D=growD  #subset(growD,allEdge==0)
  D$logarea=log(D$area)
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
    distDfile=paste("../../../speciesData/",sppList[i],"/",sppList[i],"_genet_xy.csv",sep="")
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

  library(INLA)
  #Set up ID variables for INLA random effects
  D$yearID <- D$year #for random year offset on intercept
  
  #Instead of full model, match the structure of the quadrat-based IBM regressions
  formula2 <- survives ~ logarea*crowd+
    ppt1+TmeanSpr1+ppt2+TmeanSpr2+
    f(Group, model="iid", prior="normal",param=c(1,0.001))+
    f(yearID, model="iid", prior="normal",param=c(1,0.001))+
    f(year, logarea, model="iid", prior="normal",param=c(1,0.001))
  
  outINLA <- inla(formula2, data=D,
                  family=c("binomial"), verbose=TRUE,
                  control.compute=list(dic=T,mlik=T),
                  control.predictor = list(link = 1),
                  Ntrials=rep(1,nrow(D)))
  summary(outINLA)
      
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

  # write output
  write.table(params,outfile,row.names=F,sep=",")
  
  #fixed effects
  fixedClimate <- as.data.frame(fixed[4:7,1]) 
  colnames(fixedClimate) <- "Mean"
  fixedClimate$Covariate <- rownames(fixed[4:7,])
  fixedCIs <- as.data.frame(outINLA$summary.fixed)[4:7,c(3,5)]
  fixedClimateTable <- data.frame(Covariate = fixedClimate$Covariate,
                                  Mean = fixedClimate$Mean,
                                  LowerCI = fixedCIs[,1],
                                  UpperCI = fixedCIs[,2])
  fixedEffects[[spp]] <- fixedClimateTable
  
  #Set up constant and climate models for comparison
  climate <- survives ~ logarea*crowd+
    ppt1+TmeanSpr1+ppt2+TmeanSpr2+
    f(Group, model="iid", prior="normal",param=c(1,0.001))
  
  constant <- survives ~ logarea*crowd+
    f(Group, model="iid", prior="normal",param=c(1,0.001))
  
  climateINLA <- inla(climate, data=D,
                      family=c("binomial"), verbose=TRUE,
                      control.compute=list(dic=T,mlik=T),
                      control.predictor = list(link = 1),
                      Ntrials=rep(1,nrow(D)))
  constantINLA <- inla(constant, data=D,
                       family=c("binomial"), verbose=TRUE,
                       control.compute=list(dic=T,mlik=T),
                       control.predictor = list(link = 1),
                       Ntrials=rep(1,nrow(D)))
  
  #Calculate sum of squared residuals for models
  ssrFull <- as.numeric(outINLA$dic$dic-outINLA$dic$p.ef)
  ssrConstant <- as.numeric(constantINLA$dic$dic-constantINLA$dic$p.ef)
  ssrClimate <- as.numeric(climateINLA$dic$dic-climateINLA$dic$p.ef)
  climContribution <- (ssrClimate-ssrConstant)/(ssrFull-ssrConstant)
  SSRs[spp,] <- c(length(D$logarea), ssrConstant, ssrClimate, ssrFull, climContribution) 
  
}

#Produce LATEX and HTML table of SSRs and climate covariate contributions to vital rate
ssrTab <- data.frame(Species = sppList,
                     N = SSRs[,1],
                     SSR1 = SSRs[,2],
                     SSR2 = SSRs[,3],
                     SSR3 = SSRs[,4],
                     SSRc = SSRs[,5])
colnames(ssrTab) <- c("Species", "N", "Constant model", "Climate model", 
                      "Full model", "Contribution of climate covariates")
library(xtable)
xtable(ssrTab)
print(xtable(ssrTab), type="html", file="ClimateCovariatesContributionTable.html")
write.csv(ssrTab, "ClimCovariateContribution.csv")
                                        
#Produce LATEX and HTML table of fixed climate effects means with 95% CIs
spp4tab <- (rep(as.character(sppList), each = 4))
dFixed <- as.data.frame(rbind(fixedEffects[[1]], 
                              fixedEffects[[2]],
                              fixedEffects[[3]],
                              fixedEffects[[4]]))
dFixed$Species <- spp4tab
xtable(dFixed[c(5,1,2,3,4)])
print(xtable(dFixed[c(5,1,2,3,4)]), type="html", file="fixedClimateEffectsTable.html")
print(xtable(dFixed[c(5,1,2,3,4)]), type="latex", file="fixedClimateEffectsTable.tex")

