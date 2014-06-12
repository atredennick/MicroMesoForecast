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
              clim=clim,
              nObs=length(yCov))

##CHOOSE THE MODEL
# modelFile <- "IdahoCoverIBM_JAGS.R"
# modelFile <- "IdahoCoverIBM_NoCLimate_JAGS.R"
modelFile <- "IdahoCoverIBM_JAGS_NoRandEffects.R"

jm <- jags.model(modelFile,
                data=dataJ, n.chains=1, n.adapt = 1000)
update(jm, n.iter=1000)
zm <- coda.samples(jm, variable.names=c("betaMean"), n.iter=1000, n.thin=1)

zmStat <- summary(zm)$stat
zmStat
zmQuant <- summary(zm)$quantile

# #Plot the data and model time series
# #first get means over quads
# quadAvgD <- ddply(allD, .(year), summarise,
#                   year = mean(climYear)-1,
#                   coverAvg = mean(cover, na.rm=TRUE),
#                   coverSD = sd(cover, na.rm=TRUE))
# quadAvgD$Pred <- zmStat[,1]
# quadAvgD$Low <- zmQuant[,1]
# quadAvgD$High <- zmQuant[,5]
# 
# ggplot(quadAvgD)+
#   geom_ribbon(aes(x=year, ymin=Low, ymax=High, alpha=0.5), fill="steelblue", color=NA)+
#   geom_errorbar(aes(x=year, ymax=(coverAvg+coverSD), ymin=(coverAvg-coverSD)), color="grey65")+
#   geom_point(aes(x=year, y=coverAvg), size=4)+
#   geom_line(aes(x=year, y=Pred), color="white", size=1)+
#   theme_few()+
#   ylab("Cover (%)") + xlab("Year")+
# #   scale_y_continuous(limits=c(0,12))+
#   guides(alpha=FALSE)
# 
