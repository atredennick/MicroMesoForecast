##  This is a script to be sourced from within the vital rate
##    regression fitting scripts to estimate crowding for a 
##    given dataset (e.g., one with one year left out).

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com
##  Date:   4.02.2015

##  CURRENTLY WORKING VERSION

# These are alphas, see Chu and Adler 2015 for description
alpha.effect=c(0.010,0.025,0.020,0.036)  ##take this from multiple-species IPM


####
#### Calculate crowding 
####
for(s in 1:length(sppList)){
  doSpp<-sppList[s]
  growDfile=paste("../../../speciesData/",doSpp,"/growDnoNA.csv",sep="")
  D=read.csv(growDfile)
  D$quad=as.character(D$quad)
  
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
  distD <- subset(distD, year!=year_to_leave_out)
  W=matrix(NA,dim(D)[1],length(sppList))
  crowd <- matrix(NA,dim(D)[1],length(sppList))
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
  crowd[,s] <- W[,which(sppList==doSpp)]
}#end species loop



