##  Explore PASM recruitment model
library(ggplot2)
library(plyr)
library(reshape2)

####
####  Read in recruitment parameters
####
long <- readRDS("recruitment_stanmcmc_PASM.RDS")
intercepts <- long[grep("a", long$Parameter), ] #get intercept terms
rmtheta <- grep("theta", intercepts$Parameter) #id theta rows
intercepts <- intercepts[-rmtheta, ] #remove theta rows
mu_ids <- grep("a_mu", intercepts$Parameter) #id mean intercepts
mu_intercept <- intercepts[mu_ids,] #subset mean intercepts
yr_intercept <- intercepts[-mu_ids,] #remove mean intercepts
dd <- long[grep("dd", long$Parameter),]

pdf("/Users/atredenn/Desktop/PASMrec.PDF")

####
####  Plot correlation between intercept and density-dependence
####
##  Mean intercepts
rho <- round(cor(dd$value, mu_intercept$value),2)
plot(dd$value, mu_intercept$value, main = rho)

##  Yearly intercepts
yrs <- 13
par(mfrow = c(4,4))
for(i in 1:yrs){
  param <- paste("a[", as.character(i), "]", sep="")
  tmpid <- which(yr_intercept$Parameter == param)
  tmp <- yr_intercept[tmpid, ]
  rhotmp <- round(cor(dd$value, tmp$value),2)
  plot(dd$value, tmp$value, main=rhotmp)
}


####
####  Look at raw data for all species to see if PASM is outlier
####
sppList=sort(c("BOGR","HECO","PASM","POSE"))
for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/edited/recArea.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/recArea.csv", sep=""))
    sppD$species <- doSpp 
  }
  
  sppD$Group <- as.factor(substr(sppD$quad,1,1)) #add by Chengjin
  sppD <- sppD[,c("quad","year","NRquad","totParea","Group")]
  names(sppD)[3] <- paste("R.",sppList[spp],sep="")
  names(sppD)[4] <- paste("cov.",sppList[spp],sep="")
  if(spp==1){
    D <- sppD
  }else{
    D <- merge(D,sppD,all=T)
  }
}
D[is.na(D)]=0  # replace missing values 

# calculate mean cover by group and year
tmpD <- D[,c("quad","year","Group",paste("cov.",sppList,sep=""))]
tmpD <- aggregate(tmpD[,4:NCOL(tmpD)],
                  by=list("year"=tmpD$year,"Group"=tmpD$Group),FUN=mean)
names(tmpD)[3:NCOL(tmpD)] <- paste("Gcov.",sppList,sep="")
D <- merge(D,tmpD,all.x=T)

parents1 <- as.matrix(D[,c(paste("cov.",sppList,sep=""))])/100 ##convert from absolute cover to [1,100] range
parents2 <- as.matrix(D[,c(paste("Gcov.",sppList,sep=""))])/100

##for species 1
tmp1L <- which(parents1[,1]==0) ##lcoal
tmp1G <- which(parents2[,1]==0) ##Group
tmp1 <- intersect(tmp1L,tmp1G)
##for species 2
tmp2L <- which(parents1[,2]==0)
tmp2G <- which(parents2[,2]==0)
tmp2 <- intersect(tmp2L,tmp2G)
##for species 3
tmp3L <- which(parents1[,3]==0)
tmp3G <- which(parents2[,3]==0)
tmp3 <- intersect(tmp3L,tmp3G)
##for species 4
tmp4L <- which(parents1[,4]==0)
tmp4G <- which(parents2[,4]==0)
tmp4 <- intersect(tmp4L,tmp4G)

tmp <- unique(c(tmp1,tmp2,tmp3,tmp4))

if(length(tmp)>0){
  parents1 <- parents1[-tmp,] ##remove them
  parents2 <-parents2[-tmp,] ##remove them
  y <- as.matrix(D[,c(paste("R.",sppList,sep=""))])[-tmp,] ##remove them  
  year <- D$year[-tmp] ##remove them
  Nyrs <- length(unique(D$year))
  N <- dim(D)[1]-length(tmp) ##reduce
  Nspp <- length(sppList)
  Group <- as.numeric(as.factor(D$Group))[-tmp] ##remove them ##first turn it as FACTOR, then to NUMERIC
  Ngroups <- length(unique(Group))
} else {
  y <- as.matrix(D[,c(paste("R.",sppList,sep=""))])
  year <- D$year
  Nyrs <- length(unique(D$year))
  N <- dim(D)[1]
  Nspp <- length(sppList)
  Group <- as.numeric(as.factor(D$Group)) ##first turn it as FACTOR, then to NUMERIC
  Ngroups <- length(unique(Group))
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

print(ggplot(allD, aes(parents1, recruits))+
  geom_point()+
  facet_wrap("species", scales="free")+
  ggtitle("parents1"))

print(ggplot(allD, aes(parents2, recruits))+
  geom_point()+
  facet_wrap("species", scales="free")+
  ggtitle("parents2"))

dev.off()
