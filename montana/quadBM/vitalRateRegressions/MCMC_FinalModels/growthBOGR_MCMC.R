#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

doSpp <- "BOGR"

#load libraries
library(rjags)
library(coda)
library(ggmcmc)

#bring in data
allD <- read.csv("../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
sppList <- as.character(unique(allD$Species))
sppD <- subset(allD, Species==doSpp)

climD <- read.csv("../../../weather/Climate.csv")

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


####
#### Set up data structure for JAGS
####
nGrp <- length(unique(growD$group))
nYrs <- length(unique(growD$year))
nObs <- nrow(growD)
size <- growD$percLagCover
Y <- growD$percCover
yrs <- (growD$year)-32
grp <- as.numeric(as.factor(growD$group))

####
#### Run MCMC
####
iterations <- 20000
dataJ <- list(nGrp=nGrp, nYrs=nYrs, nObs=nObs, size=size, Y=Y, yrs=yrs, grp=grp)
mod <- jags.model(paste("growth", doSpp, "_JAGS.R", sep=""), data=dataJ, n.chains=3, n.adapt=5000)
update(mod, n.iter = 10000)
out <- coda.samples(mod, c("beta", "intYr", "intG"),
                    n.iter=iterations, n.thin=10)

####
#### Check for convergence
####
gelman.diag(out)
# heidel.diag(out)
# gelman.plot(out)
# plot(out)

####
#### Convert to dataframe for export and get other summaries
####
outC <- rbind(out[[1]][(iterations-999):iterations,], 
              out[[2]][(iterations-999):iterations,], 
              out[[3]][(iterations-999):iterations,])

saveRDS(outC, file = paste(doSpp, "_GrowthParamsMCMC.rds", sep=""))
?gelman.diag()


