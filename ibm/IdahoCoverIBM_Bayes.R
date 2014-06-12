# IBM version: each "individual" is a quadrat, 
# we track cover in each quadrat over time

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

coverD <- subset(allD,!is.na(cover))  
timeN <- coverD$climYear-1927
quadN <- coverD$quadNum
nProc <- allD$climYear-1927
quadN <- allD$quadNum
y <- coverD$cover
numQuadsUnique <- length(unique(allD$quads))

dataJ <- list(timeN=timeN,
              quadN=quadN,
              nProc=length(unique(nProc)),
              nQuad=length(unique(allD$quad)),
              y=y,
              clim=clim,
              nObs=length(y))

##CHOOSE THE MODEL
# modelFile <- "IdahoCoverIBM_JAGS.R"
modelFile <- "IdahoCoverIBM_NoCLimate_JAGS.R"

jm <- jags.model(modelFile,
                data=dataJ, n.chains=1, n.adapt = 1000)
update(jm, n.iter=1000)
zm <- coda.samples(jm, variable.names=c("muN"), n.iter=1000, n.thin=1)

zmStat <- summary(zm,)$stat
zmQuant <- summary(zm)$quantile

#Plot the data and model time series
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

