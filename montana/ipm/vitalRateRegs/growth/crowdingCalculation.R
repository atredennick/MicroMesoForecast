####
#### Script to calculate crowding index for ecah species
####
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

sppList=sort(c("BOGR","HECO","PASM","POSE"))
alpha.effect=c(0.010,0.025,0.020,0.036)  ##take this from multiple-species IPM


####
#### calculate crowding 
####
for(s in 1:length(sppList)){
  doSpp<-sppList[s]
  if(doSpp == "BOGR"){
    growDfile=paste("../../../speciesData/",doSpp,"/edited/growDnoNA.csv",sep="")
  } else{
    growDfile=paste("../../../speciesData/",doSpp,"/growDnoNA.csv",sep="")
  }
  D=read.csv(growDfile)
  D$quad=as.character(D$quad)
  
  if(doSpp == "HECO"){
    ##then we moved some specific points:
    tmp2<-which(D$quad=="A12" & D$year==44)
    tmp3<-which(D$quad=="B1"  & D$year==44)
    tmp41<-which(D$quad=="E4" & D$year==33) 
    tmp42<-which(D$quad=="E4" & D$year==34) 
    tmp43<-which(D$quad=="E4" & D$year==43)
    tmp44<-which(D$quad=="E4" & D$year==44)
    
    tmpONE<-c(tmp2,tmp3,tmp41,tmp42,tmp43,tmp44)
    if(length(tmpONE)>0) D<-D[-tmpONE,]
  }
  
  for(i in 1:length(sppList)){
    if(sppList[i]=="BOGR"){
      distDfile=paste("../../../speciesData/",sppList[i],"/edited/",sppList[i],"_genet_xy_edited.csv",sep="")
    }else{
      distDfile=paste("../../../speciesData/",sppList[i],"/",sppList[i],"_genet_xy.csv",sep="")
    }
    
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
  crowdSame <- as.data.frame(W[,which(sppList==doSpp)])
  crowdSame$xID <- D$X
  colnames(crowdSame) <- c("W","X")
  write.csv(crowdSame, paste(doSpp,"growthCrowding.csv",sep=""))
}#end species loop



