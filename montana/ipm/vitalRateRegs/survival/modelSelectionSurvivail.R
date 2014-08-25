#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

sppList=sort(c("BOGR","HECO","PASM","POSE"))
alpha.effect=c(0.024,0.047,0.045,0.047) # for spp in alphabetical order

for(spp in 1:length(sppList)){
  doSpp=sppList[spp]
  climD=read.csv("../../../weather/Climate.csv")
  climD[2:ncol(climD)] <- scale(climD[2:ncol(climD)], center = TRUE, scale = TRUE)
  
  outfile=paste("Growth_params_",doSpp,".csv",sep="")
  
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
#   
#   #First do random effects structure selection by DIC
#   m1 <- survives ~ logarea*crowd+
#     ppt1*TmeanSpr1+ppt2*TmeanSpr2+
#     f(Group, model="iid", prior="normal",param=c(1,0.001))+
#     f(yearID, model="iid", prior="normal",param=c(1,0.001))+
#     f(year, logarea, model="iid", prior="normal",param=c(1,0.001))
#   
#   out1 <- inla(m1, data=D,
#                   family=c("binomial"), verbose=FALSE,
#                   control.compute=list(dic=T,mlik=T),
#                   control.predictor = list(link = 1),
#                Ntrials=rep(1,nrow(D)))
#   
#   m2 <- survives ~ logarea*crowd+
#     ppt1*TmeanSpr1+ppt2*TmeanSpr2+
#     f(Group, model="iid", prior="normal",param=c(1,0.001))+
#     f(yearID, model="iid", prior="normal",param=c(1,0.001))
#   
#   out2 <- inla(m2, data=D,
#                   family=c("binomial"), verbose=FALSE,
#                   control.compute=list(dic=T,mlik=T),
#                   control.predictor = list(link = 1),
#                Ntrials=rep(1,nrow(D)))
#                
#   m3 <- survives ~ logarea*crowd+
#     ppt1*TmeanSpr1+ppt2*TmeanSpr2+
#     f(Group, model="iid", prior="normal",param=c(1,0.001))+
#     f(year, logarea, model="iid", prior="normal",param=c(1,0.001))
#   
#   out3 <- inla(m3, data=D,
#                   family=c("binomial"), verbose=FALSE,
#                   control.compute=list(dic=T,mlik=T),
#                   control.predictor = list(link = 1),
#                Ntrials=rep(1,nrow(D)))
#   
#   m4 <- survives ~ logarea*crowd+
#     ppt1*TmeanSpr1+ppt2*TmeanSpr2+
#     f(yearID, model="iid", prior="normal",param=c(1,0.001))+
#     f(year, logarea, model="iid", prior="normal",param=c(1,0.001))
#   
#   out4 <- inla(m4, data=D,
#                   family=c("binomial"), verbose=FALSE,
#                   control.compute=list(dic=T,mlik=T),
#                   control.predictor = list(link = 1),
#                Ntrials=rep(1,nrow(D)))
#   
#   m5 <- survives ~ logarea*crowd+
#     ppt1*TmeanSpr1+ppt2*TmeanSpr2+
#     f(Group, model="iid", prior="normal",param=c(1,0.001))
#   
#   out5 <- inla(m5, data=D,
#                   family=c("binomial"), verbose=FALSE,
#                   control.compute=list(dic=T,mlik=T),
#                   control.predictor = list(link = 1),
#                Ntrials=rep(1,nrow(D)))
#   
#   m6 <- survives ~ logarea*crowd+
#     ppt1*TmeanSpr1+ppt2*TmeanSpr2+
#     f(yearID, model="iid", prior="normal",param=c(1,0.001))
#   
#   out6 <- inla(m6, data=D,
#                   family=c("binomial"), verbose=FALSE,
#                   control.compute=list(dic=T,mlik=T),
#                   control.predictor = list(link = 1),
#                Ntrials=rep(1,nrow(D)))
#   
#   m7 <- survives ~ logarea*crowd+
#     ppt1*TmeanSpr1+ppt2*TmeanSpr2+
#     f(year, logarea, model="iid", prior="normal",param=c(1,0.001))
#   
#   out7 <- inla(m7, data=D,
#                   family=c("binomial"), verbose=FALSE,
#                   control.compute=list(dic=T,mlik=T),
#                   control.predictor = list(link = 1),
#                Ntrials=rep(1,nrow(D)))
#   out1$dic$dic
#   out2$dic$dic
#   out3$dic$dic
#   out4$dic$dic
#   out5$dic$dic
#   out6$dic$dic
#   out7$dic$dic
  
  #By DIC, the first model (full mixed effects) should be used
  #So with that random effects stucture now let's do the climate effects
m1 <- survives ~ logarea*crowd+
  ppt1*TmeanSpr1+ppt2*TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out1 <- inla(m1, data=D,
             family=c("binomial"), verbose=TRUE,
             control.compute=list(dic=T,mlik=T),
             control.predictor = list(link = 1),
             Ntrials=rep(1,nrow(D)))

m2 <- survives ~ logarea*crowd+
  ppt1*TmeanSpr1+ppt2+TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out2 <- inla(m2, data=D,
             family=c("binomial"), verbose=TRUE,
             control.compute=list(dic=T,mlik=T),
             control.predictor = list(link = 1),
             Ntrials=rep(1,nrow(D)))

m3 <- survives ~ logarea*crowd+
  ppt1*TmeanSpr1+ppt2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out3 <- inla(m3, data=D,
             family=c("binomial"), verbose=TRUE,
             control.compute=list(dic=T,mlik=T),
             control.predictor = list(link = 1),
             Ntrials=rep(1,nrow(D)))

m4 <- survives ~ logarea*crowd+
  ppt1*TmeanSpr1+TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out4 <- inla(m4, data=D,
             family=c("binomial"), verbose=TRUE,
             control.compute=list(dic=T,mlik=T),
             control.predictor = list(link = 1),
             Ntrials=rep(1,nrow(D)))

m5 <- survives ~ logarea*crowd+
  ppt1*TmeanSpr1+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out5 <- inla(m5, data=D,
             family=c("binomial"), verbose=TRUE,
             control.compute=list(dic=T,mlik=T),
             control.predictor = list(link = 1),
             Ntrials=rep(1,nrow(D)))

m6 <- survives ~ logarea*crowd+
  ppt1+TmeanSpr1+ppt2*TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out6 <- inla(m6, data=D,
             family=c("binomial"), verbose=TRUE,
             control.compute=list(dic=T,mlik=T),
             control.predictor = list(link = 1),
             Ntrials=rep(1,nrow(D)))

m7 <- survives ~ logarea*crowd+
  ppt1+ppt2*TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out7 <- inla(m7, data=D,
             family=c("binomial"), verbose=TRUE,
             control.compute=list(dic=T,mlik=T),
             control.predictor = list(link = 1),
             Ntrials=rep(1,nrow(D)))

m8 <- survives ~ logarea*crowd+
  TmeanSpr1+ppt2*TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out8 <- inla(m8, data=D,
             family=c("binomial"), verbose=TRUE,
             control.compute=list(dic=T,mlik=T),
             control.predictor = list(link = 1),
             Ntrials=rep(1,nrow(D)))

m9 <- survives ~ logarea*crowd+
  ppt2*TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out9 <- inla(m9, data=D,
             family=c("binomial"), verbose=TRUE,
             control.compute=list(dic=T,mlik=T),
             control.predictor = list(link = 1),
             Ntrials=rep(1,nrow(D)))

m10 <- survives ~ logarea*crowd+
  ppt1+TmeanSpr1+ppt2+TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out10 <- inla(m10, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m11 <- survives ~ logarea*crowd+
  ppt1+TmeanSpr1+ppt2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out11 <- inla(m11, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m12 <- survives ~ logarea*crowd+
  ppt1+TmeanSpr1+TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out12 <- inla(m12, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m13 <- survives ~ logarea*crowd+
  ppt1+TmeanSpr1+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out13 <- inla(m13, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m14 <- survives ~ logarea*crowd+
  ppt1+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out14 <- inla(m14, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m15 <- survives ~ logarea*crowd+
  TmeanSpr1+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out15 <- inla(m15, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m16 <- survives ~ logarea*crowd+
  ppt1+ppt2+TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out16 <- inla(m16, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m17 <- survives ~ logarea*crowd+
  TmeanSpr1+ppt2+TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out17 <- inla(m17, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m18 <- survives ~ logarea*crowd+
  ppt2+TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out18 <- inla(m18, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m19 <- survives ~ logarea*crowd+
  ppt2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out19 <- inla(m19, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m20 <- survives ~ logarea*crowd+
  TmeanSpr2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out20 <- inla(m20, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m21 <- survives ~ logarea*crowd+
  TmeanSpr2+ppt1+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out21 <- inla(m21, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m22 <- survives ~ logarea*crowd+
  TmeanSpr1+ppt2+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out22 <- inla(m22, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m23 <- survives ~ logarea*crowd+
  TmeanSpr2+TmeanSpr1+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out23 <- inla(m23, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))
m24 <- survives ~ logarea*crowd+
  ppt2+ppt1+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out24 <- inla(m24, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

m25 <- survives ~ logarea*crowd+
  f(Group, model="iid", prior="normal",param=c(1,0.001))+
  f(yearID, model="iid", prior="normal",param=c(1,0.001))+
  f(year, logarea, model="iid", prior="normal",param=c(1,0.001))

out25 <- inla(m25, data=D,
              family=c("binomial"), verbose=TRUE,
              control.compute=list(dic=T,mlik=T),
              control.predictor = list(link = 1),
              Ntrials=rep(1,nrow(D)))

dics <- c(out1$dic$dic, out2$dic$dic, out3$dic$dic, out4$dic$dic, out5$dic$dic,
          out6$dic$dic, out7$dic$dic, out8$dic$dic, out9$dic$dic, out10$dic$dic,
          out11$dic$dic, out12$dic$dic, out13$dic$dic, out14$dic$dic, out15$dic$dic,
          out16$dic$dic, out17$dic$dic, out18$dic$dic, out19$dic$dic, out20$dic$dic,
          out21$dic$dic, out22$dic$dic, out23$dic$dic, out24$dic$dic, out25$dic$dic)
dics <- as.data.frame(dics)
dics$model <- c(1:25)
dicRank <- dics[with(dics, order(dics)), ] #Can make this into table (top 5?) post hoc
finModTxt <- paste("out", as.character(dicRank[1,2]), sep="")
finMod <- eval(parse(text=finModTxt)) 
summary(finMod)

write.csv(dicRank, file = paste(sppList[spp], "_DICranks.csv"))
save(finMod, file = paste(sppList[spp], "_finalSurvivalModel.Rdata"))
}
