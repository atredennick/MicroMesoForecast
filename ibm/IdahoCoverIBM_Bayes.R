# IBM version: each "individual" is a quadrat, 
# we track cover in each quadrat over time

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#


#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(reshape2)
library(plyr)
library(ggplot2)
library(ggthemes)
library(rjags)
library(coda)

doSpp="Artemisia tripartita"
# doSpp="Poa secunda"
# doSpp="Pseudoroegneria spicata"

# 1. Import and format data ------------------------------------
climD=read.csv("Climate_full.csv")
climD=subset(climD,year>1927)  # remove early years with NAs

# sum focal species cover by quad and year, make sure to get zeros
rawD=read.csv("allrecords_cover.csv")
QuadYrList=unique(rawD[,c("quad","year")],MARGIN=2)
rawD=subset(rawD, species==doSpp)
quadD=aggregate(area~quad+year,data=rawD,FUN=sum)
quadD$cover=quadD$area*100
quadD=merge(quadD,QuadYrList,all.y=T) 
quadD$cover[is.na(quadD$cover)]=0  # add in zeros

# create lag cover  variable
tmp=quadD[,c("quad","year","cover")]
tmp$year=tmp$year+1
names(tmp)[3]="lag.cover"
quadD=merge(quadD,tmp,all.x=T)

#add in missing data so that each quadrat is represented each climate year, even if cover is NA
allYears <- seq(min(unique(quadD$year)), max(unique(quadD$year)), 1)
allQuads <- unique(quadD$quad)
allQuadYears <- data.frame(quad = sort(rep(allQuads, length(allYears))),
                           year = rep(allYears, length(allQuads)))
allD <- merge(allQuadYears, quadD, by=c("quad", "year"), all.x=T)

# merge in climate data
allD$climYear=allD$year+1900-1  
allD=merge(allD,climD,by.x="climYear",by.y="year")
allD$quadNum <- as.numeric(substr(allD$quad,2,3))


# 2. Set up data structures and indexing for JAGS model ----------------
# climDs=subset(climD,year<max(allD$climYear))
clim <- allD[7:(ncol(x=allD)-1)]

### FOR LIKELIHOODS ###
#All cover observations
coverD <- subset(allD,!is.na(cover))  
timeNcov <- coverD$climYear-1927
quadNcov <- coverD$quadNum
yCov <- coverD$cover

#Colonization observations
colD <- subset(coverD,lag.cover==0)
yCol <- ifelse(colD$cover>0,1,0)
timeNcol <- colD$climYear-1927
quadNcol <- colD$quadNum

#Survival observations
survD <- subset(coverD,lag.cover>0)
ySurv <- ifelse(survD$cover>0,1,0)
timeNsurv <- survD$climYear-1927
quadNsurv <- survD$quadNum

#Growth observations
growD <- subset(coverD,lag.cover>0 & cover>0)
yGrow <- growD$cover
timeNgrow <- growD$climYear-1927
quadNgrow <- growD$quadNum

### FOR PROCESS MODEL ###
nProc <- allD$climYear-1927
quadN <- allD$quadNum
numQuadsUnique <- length(unique(allD$quads))

dataJ <- list(timeN=timeNcov,
              quadN=quadNcov,
              nProc=length(unique(nProc)),
              nQuad=length(unique(allD$quad)),
              y=yCov,
              yCol=yCol,
              timeNcol=timeNcol,
              quadNcol=quadNcol,
              ySurv=ySurv,
              timeNsurv=timeNsurv,
              quadNsurv=quadNsurv,
              yGrow=yGrow,
              timeNgrow=timeNgrow,
              quadNgrow=quadNgrow,
              clim=clim,
              nObs=length(yCov),
              nCol=length(yCol),
              nSurv=length(ySurv),
              nGrow=length(yGrow))

##CHOOSE THE MODEL
# modelFile <- "IdahoCoverIBM_JAGS.R"
# modelFile <- "IdahoCoverIBM_NoCLimate_JAGS.R"
modelFile <- "IdahoCoverIBM_JAGS_NoRandEffects.R"

n.Adapt <- 1000
n.Up <- 1000
n.Samp <- 2000

jm <- jags.model(modelFile,
                data=dataJ, n.chains=1, n.adapt = n.Adapt)
update(jm, n.iter=n.Up)
# Get sums of square residuals
zm <- coda.samples(jm, variable.names=c("fit", "bpvalue"), n.iter=n.Samp, n.thin=10)
ssq <- summary(zm)$stat[2,1]
bP <- summary(zm)$stat[1,1]

#Get coefficients
zm <- coda.samples(jm, variable.names=c("betaClim", "alphaClim", "gammaClim"), n.iter=n.Samp, n.thin=10)

zmD <- as.data.frame(zm[[1]])
zmM <- melt(zmD)
zmM$Vital <- c(rep("Survival",n.Adapt*4), rep("Growth",n.Adapt*4), rep("Colonization",n.Adapt*4))
zmM$Covariate <- rep(c(rep("ppt2",n.Adapt), rep("ppt1",n.Adapt), rep("TmeanSpr2",n.Adapt), rep("TmeanSpr1",n.Adapt)), 3)

zmStat <- summary(zm)$stat
zmStat
zmQuant <- summary(zm)$quantile

#Plot climate coefficients by vital rate
climCoefD <- data.frame(Vital = c(rep("Survival",4), rep("Growth",4), rep("Colonization",4)),
                        Covariate = rep(c("ppt2", "ppt1", "TmeanSpr2", "TmeanSpr1"), 3),
                        Mean = zmStat[,1],
                        Up = zmQuant[,5],
                        Down = zmQuant[,1],
                        Midup = zmQuant[,4],
                        Middown = zmQuant[,2])
ggplot()+
  geom_hline(yintercept=0, alpha=0.5, linetype="dashed")+
  geom_point(data=zmM,aes(x=Covariate, y=value, color=Vital), shape="-", alpha=0.01, size=4)+
  geom_point(data=climCoefD, aes(x=Covariate, y=Mean, color=Vital), fill="white", size=4, shape=21)+
#   geom_point(data=climCoefD, aes(x=Covariate, y=Mean), size=3, shape=1)+
  facet_grid(.~Vital)+
  guides(color=FALSE)+
  theme_few()

ggplot(data=climCoefD)+
  geom_errorbar(aes(x=Covariate, ymin=Down, ymax=Up, color=Vital), width=0.1)+
  geom_errorbar(aes(x=Covariate, ymin=Middown, ymax=Midup, color=Vital), width=0.000001, size=2)+
  geom_point(aes(x=Covariate, y=Mean, color=Vital), size=5)+
  geom_point(aes(x=Covariate, y=Mean), size=5, shape=1)+
  geom_hline(yintercept=0, alpha=0.5, linetype="dashed")+
  facet_grid(.~Vital)+
  guides(color=FALSE)+
  theme_few()



#Plot the data and model time series
zm <- coda.samples(jm, variable.names=c("muN"), n.iter=n.Samp, n.thin=10)
zmQuant <- summary(zm)$quantile
zmStat <- summary(zm)$stat

#first get means over quads
quadAvgD <- ddply(allD, .(year), summarise,
                  year = mean(climYear)-1,
                  coverAvg = mean(cover, na.rm=TRUE),
                  coverSD = sd(cover, na.rm=TRUE))
quadAvgD$Pred <- zmStat[,1]
quadAvgD$Low <- zmQuant[,1]
quadAvgD$High <- zmQuant[,5]

ggplot(quadAvgD)+
  geom_ribbon(aes(x=year, ymin=Low, ymax=High, alpha=0.5), fill="steelblue", color=NA)+
  geom_errorbar(aes(x=year, ymax=(coverAvg+coverSD), ymin=(coverAvg-coverSD)), color="grey65")+
  geom_point(aes(x=year, y=coverAvg), size=4)+
  geom_line(aes(x=year, y=Pred), color="white", size=1)+
  theme_few()+
  ylab("Cover (%)") + xlab("Year")+
#   scale_y_continuous(limits=c(0,12))+
  guides(alpha=FALSE)


## Plot the climate record
recClimD <- subset(x=climD, subset=year<max(allD$climYear))
fullClimD <- melt(recClimD, id.vars="year")
ggplot(data=fullClimD, aes(x=year, y=value))+
  geom_line(aes(color=variable))+
  geom_vline(xintercept=1945)+
  facet_grid(variable~., scales = "free")+
  guides(color=FALSE)+
  theme_few()


