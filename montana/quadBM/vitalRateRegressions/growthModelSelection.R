# IBM version: each "individual" is a quadrat, 
# we track cover in each quadrat over time

#This script is for model selection based on DIC, regardless of 95% CIs

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

####
#### POSSIBLE FIXED-EFFECTS -----------------------
####
source("fixedEffectsModels_Growth.R")


####
#### PRELIMINARIES --------------------------------
####

#load libraries
library(INLA)

#bring in data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../weather/Climate.csv")

####
#### LOOP SPECIES ---------------------------------
####
dics <- matrix(ncol=length(sppList), nrow=length(fixed.forms))
for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  sppD <- subset(allD, Species==doSpp)
  
  # create lag cover  variable
  tmp=sppD[,c("quad","year","totCover")]
  tmp$year=tmp$year+1
  names(tmp)[3]="lag.cover"
  sppD=merge(sppD,tmp,all.x=T)
  
  # merge in climate data
  sppD$climYear=sppD$year+1900-1  
  sppD=merge(sppD,climD,by.x="climYear",by.y="year")
  
  #Growth observations
  growD <- subset(sppD,lag.cover>0 & totCover>0)
  growD$yearID <- growD$year #for random year offset on intercept
  growD$group <- substring(growD$quad, 1, 1)
  growD$percCover <- growD$totCover/10000
  growD$percLagCover <- growD$lag.cover/10000

  #set up unchanging random effects
  random <- "+f(yearID, model='iid', prior='normal',param=c(0,0.001))+
              f(year, percLagCover, model='iid', prior='normal', param=c(0,0.001))+
              f(group, model='iid', prior='normal',param=c(0,0.001))"
  
  #Loop through climate effects
  for(j in 1:length(fixed.forms)){
    fixed <- fixed.forms[[j]]
    formula1 <- as.formula(paste(fixed, random)) 
    out1 <- inla(formula1, data=growD,
                 family=c("beta"), verbose=TRUE,
                 control.compute=list(dic=T,mlik=T),
                 control.inla = list(h = 1e-10))
    dics[j,spp] <- out1$dic$dic
  }#end model loop
}#end species loop
  

####
#### EXPORT DICs --------------------------------
####
colnames(dics) <- sppList
dics <- as.data.frame(dics)
dics$modelNum <- c(1:length(fixed.forms))
library(reshape2)
dicM <- melt(dics, id.vars = "modelNum")
colnames(dicM)[2:3] <- c("Species", "DIC")

outfile <- "./growthDIC.csv"
write.csv(dicM, outfile)

 