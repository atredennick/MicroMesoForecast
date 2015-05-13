#Quad-Based Model simulations

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

library(reshape2)
library(plyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(EnvStats)
library(rstan)
library(ggmcmc)

#bring in data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
head(scale(allD$percCover, center=TRUE, scale=TRUE))
sppList <- as.character(unique(allD$Species))

#bring in climate data
climD <- read.csv("../../weather/Climate.csv")
climD[2:6] <- scale(climD[2:6], center = TRUE, scale = TRUE)

#perturb climate data
# climD <- read.csv("../../weather/Climate.csv")
# climScale <- scale(climD[3:6], center = TRUE, scale = TRUE)
# climAvg <- apply(X = climD, MARGIN = 2, FUN = mean)
# climSD <- apply(X = climD, MARGIN = 2, FUN = sd)
# climD[c(3,5)] <- climD[c(3,5)]+(climD[c(3,5)]*0.01)
# climD[c(4,6)] <- climD[c(4,6)]+(climD[c(4,6)]*0.01)
# climD[3] <- (climD[3] - climAvg[3])/climSD[3]
# climD[4] <- (climD[4] - climAvg[4])/climSD[4]
# climD[5] <- (climD[5] - climAvg[5])/climSD[5]
# climD[6] <- (climD[6] - climAvg[6])/climSD[6]

##  Load vital rate parameters
do_species <- "POSE"
fitlong <- readRDS(paste("../vitalRateRegressions/truncNormModel/mcmc_popgrowth_", 
                         do_species, ".RDS", sep=""))
fitlong$keep <- "no"
keepseq <- seq(from = 1, to = nrow(fitlong), by = 10)
fitlong[keepseq,"keep"] <- "yes"
fitthin <- subset(fitlong, keep=="yes")

##  Break up MCMC into regression components
# Climate effects
climeff <- fitthin[grep("b2", fitthin$Parameter),]

# Yearly cover (size) effects
coveff <- fitthin[grep(glob2rx("b1[*]"), fitthin$Parameter),]
coveff$yearid <- substr(coveff$Parameter, 4, length(coveff$Parameter))
coveff$yearid <- unlist(strsplit(coveff$yearid, split=']'))

# Yearly intercepts
intercept <- fitthin[grep("a", fitthin$Parameter),]
intercept <- subset(intercept, Parameter!="a_mu")
intercept <- subset(intercept, Parameter!="tau")
intercept$yearid <- substr(intercept$Parameter, 3, length(intercept$Parameter))
intercept$yearid <- unlist(strsplit(intercept$yearid, split=']'))

# Lognormal sigma (called tau here)
tau <- fitthin[grep("tau", fitthin$Parameter),]

##  Define population growth function
growFunc <- function(N, int, slope, clims, climcovs, tau){
  mu <- int+slope*log(N)+sum(clims*climcovs)
  newN <- rlnormTrunc(1, meanlog = mu, sdlog = tau, min = 0, max = 1)
  return(newN)
}


##  Run simulations
outD <- data.frame(cover=NA, species=NA, year=NA)
tsims <- 100
cover <- numeric(tsims)
cover[1] <- 0.01
for(t in 2:tsims){
  randchain <- sample(x = climeff$Chain, size = 1)
  randiter <- sample(x = climeff$Iteration, size = 1)
  randyear <- sample(x = intercept$yearid, size = 1)
  inttmp <- subset(intercept, Chain==randchain & 
                              Iteration==randiter &
                              yearid==randyear)
  slopetmp <- subset(coveff, Chain==randchain & 
                             Iteration==randiter &
                             yearid==randyear)
  tmpclim <- subset(climeff, Chain==randchain & 
                             Iteration==randiter)
  tmptau <- subset(tau, Chain==randchain & 
                        Iteration==randiter)
  climyear <- sample(c(1:nrow(climD)), size = 1)
  climcovs <- climD[climyear,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
  climcovs$inter1 <- climcovs$ppt1*climcovs$TmeanSpr1
  climcovs$inter2 <- climcovs$ppt2*climcovs$TmeanSpr2
  cover[t] <- growFunc(N = cover[t-1], int = inttmp$value, 
                       slope = slopetmp$value, clims = tmpclim$value,
                       climcovs = climcovs, tau = tmptau$value) 
}

plot(c(1:tsims), cover*100, type="l")
