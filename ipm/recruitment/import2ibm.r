# all climate variables used
climNames=names(climD)

# import parameters
# recruitment parameters
Rpars=list(intcpt.mu=rep(0,Nspp),intcpt.yr=matrix(0,Nyrs,Nspp),betaC=matrix(0,Nspp,length(climNames)),
  intcpt.gr=matrix(0,6,Nspp),g.tau=rep(NA,Nspp),
  dd=rep(NA,Nspp),theta=rep(NA,Nspp),sizeMean=rep(NA,Nspp),sizeVar=rep(NA,Nspp),
  recSizes=list(1))

 infile="recruitment/recruit_params_Ryr.csv"    
 Rdata=read.csv(infile)
 # subset out non-essential parameters
 tmp=c(grep("lambda",row.names(Rdata)),grep("deviance",row.names(Rdata)),grep("DIC",row.names(Rdata)))   #group stuff?
 Rdata=Rdata[-tmp,]
 tmp=paste("Rpars$",row.names(Rdata),"<-",Rdata[,1],sep="")
 eval(parse(n=dim(Rdata)[1],text=tmp))
 for(i in 1:Nspp){
   infile=paste("H:/idahochart/ipm/speciesData/",sppList[i],"/recSize.csv",sep="")
   recSize=read.csv(infile)
   Rpars$sizeMean[i]=mean(log(recSize$area))
   Rpars$sizeVar[i]=var(log(recSize$area))
   #Rpars$recSizes[[i]]=recSize$area
 }
Rpars$dd=t(Rpars$dd) # c[i,j] = effect of j on i
rm(Rdata)

# recruitment
recruit=function(Rpars,sizes,spp,weather,doYear,lastID,L,expand){
    # dd is a matrix of dim Nspp X Nspp
    # sizes and spp are vectors of the same length (=N plants)
    
    # calculate total areas
    totArea=aggregate(sizes,by=list("spp"=spp),FUN=sum)
    # put in missing zeros
    tmp=data.frame("spp"=1:length(sppList))
    totArea=merge(totArea,tmp,all.y=T)
    totArea[is.na(totArea)]=0
    totArea=totArea[,2]/((L*expand)^2)*100  # scale to % cover   
    # calculate recruits
    lambda=rep(NA,Nspp) # seed production
    for(i in 1:Nspp){
       lambda[i]=totArea[i]*exp(Rpars$intcpt.yr[doYear,i]+sum(Rpars$betaC[i,]*weather)+sqrt(totArea[i])*Rpars$dd[i]) 
    }
    # number of draws from distribution depends on size of landscape
    NN=rnbinom(length(lambda)*expand^2,mu=lambda,size=Rpars$theta)  
    NN=rowSums(matrix(NN,length(lambda),expand^2))      
    x=y=spp=id=size=NULL
    for(i in 1:Nspp){
      if(NN[i]>0){
        #get recruit sizes 
        size=c(size,exp(rnorm(NN[i],Rpars$sizeMean[i],sqrt(Rpars$sizeVar[i]))))
        if(sum(is.na(size))>0) stop("Check recruit sizes")
        #assign random coordinates
        x=c(x,expand*L*runif(NN[i])); y=c(y,expand*L*runif(NN[i]))
        spp=c(spp,rep(i,NN[i]))
        #assign genet ID's
        if(length(id)>0) lastID=max(id)
        id=c(id,(1:NN[i]+lastID))
      }      
   } # next i
   # output
   size[size<minSize]=minSize
   out=cbind(spp,size,x,y,id)
   return(out)
}
