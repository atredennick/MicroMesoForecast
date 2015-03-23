##  This script runs leave-one-year-out fit for
##    the truncNorm quad-based population growth model.
##    Set 'year_index_for_leave_out' for the year to NOT use
##    for fitting. Model parameters are then used to estimate
##    cover for the unfitted year as validation. We do this
##    for all 13 observation years.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

year_index_for_leave_out <- 11

#load libraries
library(rjags)
library(coda)
load.module("dic")

#bring in data
allD <- read.csv("../../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))

climD <- read.csv("../../../../weather/Climate.csv")
climD[2:6] <- scale(climD[2:6], center = TRUE, scale = TRUE)

backD <- data.frame(climYear=NA,
                    quad = NA,
                    year= NA,
                    totCover= NA,
                    Species= NA,
                    propCover= NA,
                    lag.cover= NA,
                    pptLag= NA,
                    ppt1= NA,
                    TmeanSpr1= NA,
                    ppt2= NA,
                    TmeanSpr2= NA,
                    TmeanSum1= NA,
                    TmeanSum2= NA,
                    yearID= NA,
                    group = NA,
                    percCover = NA,
                    percLagCover = NA)

#loop through species and remake data frame
for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  sppD <- subset(allD, Species==doSpp)
  
  # create lag cover variable
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
  backD <- rbind(backD, growD)
}#end species loop

growD <- backD[2:nrow(backD),]

##  Leave one year out for fitting
all_years <- unique(growD$year)
year_to_leave_out <- all_years[year_index_for_leave_out]
growD <- subset(growD, year!=year_to_leave_out)

####
#### Set up data structure for JAGS
####
nGrp <- length(unique(growD$group))
nYrs <- length(unique(growD$year))
nObs <- nrow(growD)
C <- growD$percCover

#this is to get the sequencing right for indexing in the JAGS model
yrs <- as.numeric(as.factor((growD$year)-32))

grp <- as.numeric(as.factor(growD$group))
TmeanSpr1 <- growD$TmeanSpr1
TmeanSpr2 <- growD$TmeanSpr2
ppt1 <- growD$ppt1
ppt2 <- growD$ppt2
pptLag <- growD$pptLag
nSpp <- length(unique(growD$Species))
spp <- as.numeric(as.factor(growD$Species))
X <- log(growD$percLagCover)

####
#### Run MCMC
####
iterations <- 50000
adapt <- 10000
dataJ <- list(nGrp=nGrp, nYrs=nYrs, nObs=nObs, C=C, X=X, yrs=yrs, grp=grp,
              TmeanSpr1=TmeanSpr1, TmeanSpr2=TmeanSpr2, ppt1=ppt1, ppt2=ppt2, 
              pptLag=pptLag,spp=spp, nSpp=nSpp)
mod <- jags.model("growthAllSpp_JAGS.R", data=dataJ, n.chains=3, n.adapt=adapt)
update(mod, n.iter = (iterations))
out <- coda.samples(mod, c("intYr", "intercept", "beta", "betaSpp",
                           "intG", "temp1", "temp2", "rain1", "rain2",
                           "rainLag", "tau"),
                    n.iter=iterations, n.thin=10)


####
#### Check for convergence
####
gelmDiag <- gelman.diag(out)

####
#### Convert to dataframe for export and get other summaries
####
outC <- rbind(out[[1]][(iterations-999):iterations,], 
              out[[2]][(iterations-999):iterations,], 
              out[[3]][(iterations-999):iterations,])

outStat <- as.data.frame(summary(out)$stat)
outQuant <- as.data.frame(summary(out)$quantile)

sppNames <- c(rep(sppList, 12+6+12+4+4))
outStat$species <- sppNames
outQuant$species <- sppNames

uniq_years <- unique(growD$year)
year_names <- c(rep(uniq_years, each=4), rep(NA,7*4),
                rep(uniq_years, each=4), rep(NA,7*4))
outStat$year <- year_names
outQuant$year <- year_names

saveRDS(outC, file = paste("growthParamsMCMC_", year_to_leave_out, ".rds", sep=""))
write.csv(gelmDiag[[1]], file=paste("growthGelman_", year_to_leave_out, ".csv", sep=""))
write.csv(outStat, file=paste("growthStats_", year_to_leave_out, ".csv", sep=""))
write.csv(outQuant, file=paste("growthQuants_", year_to_leave_out, ".csv", sep=""))



