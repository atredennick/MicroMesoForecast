#Quad-Based Model simulations

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))
climvec <- readRDS("../../../climate_year_sequence.rds")
do_species <- "BOGR"
tsims <- 2500
burn.in <- 500
outs <- matrix(ncol=2, nrow=3)

library(reshape2)
library(plyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(EnvStats)

#bring in data
allD <- read.csv("../../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
# head(scale(allD$percCover, center=TRUE, scale=TRUE))
sppList <- as.character(unique(allD$Species))

#perturb climate data
climD <- read.csv("../../../weather/Climate.csv")
climScale <- scale(climD[2:6], center = TRUE, scale = TRUE)
climAvg <- apply(X = climD, MARGIN = 2, FUN = mean)
climSD <- apply(X = climD, MARGIN = 2, FUN = sd)
climD[2] <- (climD[2] - climAvg[2])/climSD[2]
climD[3] <- (climD[3] - climAvg[3])/climSD[3]
climD[4] <- (climD[4] - climAvg[4])/climSD[4]
climD[5] <- (climD[5] - climAvg[5])/climSD[5]
climD[6] <- (climD[6] - climAvg[6])/climSD[6]

fitlong <- readRDS(paste("../../vitalRateRegressions/truncNormModel/popgrowth_stanmcmc_", 
                         do_species, ".RDS", sep=""))
fitsumm <- ddply(fitlong, .(Parameter), summarise,
                 mean_value = mean(value),
                 hi_value = quantile(value, probs = 0.875),
                 lo_value = quantile(value, probs = 0.125))
fitthin <- melt(fitsumm, id.vars = "Parameter")
 
##  Break up MCMC into regression components
# Climate effects
climeff <- fitthin[grep("b2", fitthin$Parameter),]
climeff$id <- substr(climeff$Parameter, 4, length(climeff$Parameter))
climeff$id <- as.numeric(unlist(strsplit(climeff$id, split=']')))  

# Mean cover (size) effects
coveff <- fitthin[grep(glob2rx("b1_mu"), fitthin$Parameter),]
  
# Mean intercept
intercept <- fitthin[grep("a_mu", fitthin$Parameter),]
  
# Lognormal sigma (called tau here)
tau <- fitthin[grep("tau", fitthin$Parameter),]
  
##  Define population growth function
growFunc <- function(N, int, slope, clims, climcovs, tau){
  mu <- int+slope*log(N)+sum(clims[1:7]*climcovs)+sum(clims[8:12]*log(N)*climcovs[1:5])
  newN <- rlnormTrunc(1, meanlog = mu, sdlog = tau, min = 0, max = 1)
  return(newN)
}

var_names <- unique(fitthin$variable)

####
##  Simulate with mean parameters only ------------------------
####
cover <- numeric(tsims)
cover[1] <- 0.01
randset <- "mean_value"
inttmp <- subset(intercept, variable==randset)
slopetmp <- subset(coveff, variable==randset)
tmpclim <- subset(climeff, variable==randset)
tmpclim <- tmpclim[with(tmpclim, order(id)),]
tmptau <- subset(tau, variable==randset)
pb <- txtProgressBar(min=2, max=tsims, char="+", style=3, width=65)
for(t in 2:tsims){
  climyear <- climvec[t] - min(climvec)+1
  climcovs <- climD[climyear,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
  climcovs$inter1 <- climcovs$ppt1*climcovs$TmeanSpr1
  climcovs$inter2 <- climcovs$ppt2*climcovs$TmeanSpr2
  cover[t] <- growFunc(N = cover[t-1], int = inttmp$value, 
                       slope = slopetmp$value, clims = tmpclim$value,
                       climcovs = climcovs, tau = tmptau$value) 
#   if(cover[t]==0) cover[t] <- 0.01
  setTxtProgressBar(pb, t)
}#end simulation loop
plot(cover, type="l")
mean(cover)

####
##  Vary params every time step -------------------------------
####
cover <- numeric(tsims)
cover[1] <- 0.01
pb <- txtProgressBar(min=2, max=tsims, char="+", style=3, width=65)
for(t in 2:tsims){
  randchain <- sample(x = c(1:3), size = 1)
  randset <- var_names[randchain]
  inttmp <- subset(intercept, variable==randset)
  slopetmp <- subset(coveff, variable==randset)
  tmpclim <- subset(climeff, variable==randset)
  tmpclim <- tmpclim[with(tmpclim, order(id)),]
  tmptau <- subset(tau, variable==randset)
  climyear <- climvec[t] - min(climvec)+1
  climcovs <- climD[climyear,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
  climcovs$inter1 <- climcovs$ppt1*climcovs$TmeanSpr1
  climcovs$inter2 <- climcovs$ppt2*climcovs$TmeanSpr2
  cover[t] <- growFunc(N = cover[t-1], int = inttmp$value, 
                         slope = slopetmp$value, clims = tmpclim$value,
                         climcovs = climcovs, tau = tmptau$value) 
#   if(cover[t]==0) cover[t] <- 0.01
  setTxtProgressBar(pb, t)
}#end simulation loop
plot(cover, type="l")
outs[1,1] <- mean(cover[burn.in:tsims])
outs[1,2] <- sd(cover[burn.in:tsims])


####
##  Sample all params each timestep --------------------------------
####
cover <- numeric(tsims)
cover[1] <- 0.01
pb <- txtProgressBar(min=2, max=tsims, char="+", style=3, width=65)
for(t in 2:tsims){
  samps <- c(1:3)
  covtmp <- numeric(length(samps))
  for(i in samps){
    randset <- var_names[i]
    inttmp <- subset(intercept, variable==randset)
    slopetmp <- subset(coveff, variable==randset)
    tmpclim <- subset(climeff, variable==randset)
    tmpclim <- tmpclim[with(tmpclim, order(id)),]
    tmptau <- subset(tau, variable==randset)
    climyear <- climvec[t] - min(climvec)+1
    climcovs <- climD[climyear,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
    climcovs$inter1 <- climcovs$ppt1*climcovs$TmeanSpr1
    climcovs$inter2 <- climcovs$ppt2*climcovs$TmeanSpr2
    covtmp[i] <- growFunc(N = cover[t-1], int = inttmp$value, 
                         slope = slopetmp$value, clims = tmpclim$value,
                         climcovs = climcovs, tau = tmptau$value)
#     if(covtmp[i]==0) covtmp[i] <- 0.01
  }
  cover[t] <- mean(covtmp)
  setTxtProgressBar(pb, t)
}#end simulation loop
plot(cover, type="l")
outs[2,1] <- mean(cover[burn.in:tsims])
outs[2,2] <- sd(cover[burn.in:tsims])


####
##  Separate simulations for each parameter set ------------------------------
####
cover <- matrix(nrow=tsims, ncol=3)
cover[1,] <- 0.01
pb <- txtProgressBar(min=2, max=tsims, char="+", style=3, width=65)
for(i in 1:3){
  randset <- var_names[i]
  inttmp <- subset(intercept, variable==randset)
  slopetmp <- subset(coveff, variable==randset)
  tmpclim <- subset(climeff, variable==randset)
  tmpclim <- tmpclim[with(tmpclim, order(id)),]
  tmptau <- subset(tau, variable==randset)
  for(t in 2:tsims){
    climyear <- climvec[t] - min(climvec)+1
    climcovs <- climD[climyear,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
    climcovs$inter1 <- climcovs$ppt1*climcovs$TmeanSpr1
    climcovs$inter2 <- climcovs$ppt2*climcovs$TmeanSpr2
    cover[t,i] <- growFunc(N = cover[t-1,i], int = inttmp$value, 
                         slope = slopetmp$value, clims = tmpclim$value,
                         climcovs = climcovs, tau = tmptau$value) 
#     if(cover[t,i]==0) cover[t,i] <- 0.01
  }
  setTxtProgressBar(pb, t)
}
covout <- apply(cover,1, mean)
plot(covout, type="l")
outs[3,1] <- mean(covout[burn.in:tsims])
outs[3,2] <- sd(covout[burn.in:tsims])

