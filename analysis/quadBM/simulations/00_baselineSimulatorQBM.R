#Quad-Based Model simulations

#clear everything, just to be safe 
rm(list=ls(all=TRUE))


####
####  Load Libraries
####
library(reshape2);library(plyr);library(ggplot2);library(ggthemes)
library(gridExtra);library(EnvStats) # for truncated log normal functions
library(rstan);library(ggmcmc)



####
####  Preliminaries
####
do_species <- "HECO"
clim_scalers <- subset(readRDS("../../qbm_all_clim_scalers.RDS"), 
                       species==do_species & yearout==1) # get focal spp and all years (1)



####
####  Read in Data
####
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
sppList <- as.character(unique(allD$Species))

##  Bring in climate data
climD <- read.csv("../../weather/Climate.csv")

##  Load vital rate parameters
fitlong <- readRDS(paste("../vitalRateRegressions/truncNormModel/popgrowth_stanmcmc_", 
                         do_species, ".RDS", sep=""))

##  Keep every 10th MCMC iteration for random parameter draws
fitlong$keep <- "no"
keepseq <- seq(from = 1, to = nrow(fitlong), by = 10)
fitlong[keepseq,"keep"] <- "yes"
fitthin <- subset(fitlong, keep=="yes")

##  Break up MCMC into regression components
# Climate effects
climeff <- fitthin[grep("b2", fitthin$Parameter),]

# Yearly cover (size) effects
coveff <- fitthin[grep("b1[.]", fitthin$Parameter),]
coveff$yearid <- substr(coveff$Parameter, 4, length(coveff$Parameter))
coveff$yearid <- as.numeric(unlist(strsplit(coveff$yearid, split="[.]")))

# Yearly intercepts
intercept <- fitthin[grep("a", fitthin$Parameter),]
intercept <- subset(intercept, Parameter!="a_mu")
intercept <- subset(intercept, Parameter!="tau")
intercept$yearid <- substr(intercept$Parameter, 3, length(intercept$Parameter))
intercept$yearid <- as.numeric(unlist(strsplit(intercept$yearid, split='[.]')))

# Lognormal sigma (called tau here)
tau <- fitthin[grep("tau", fitthin$Parameter),]



####
####  Define Population Growth Function
####
growFunc <- function(N, int, slope, clims, climcovs, tau){
  mu <- int+slope*log(N)+sum(clims*climcovs)
  newN <- rlnormTrunc(1, meanlog = mu, sdlog = tau, min = 0, max = 1)
  return(newN)
}



####
####  Run Simulation
####
outD <- data.frame(cover=NA, species=NA, year=NA)
tsims <- 1000
cover <- numeric(tsims)
cover[1] <- 0.01 # initial cover, arbitrarily low

# Start progress bar
pb <- txtProgressBar(min=2, max=tsims, char="+", style=3, width=65)

for(t in 2:tsims){
  randchain <- sample(x = unique(climeff$Chain), size = 1) # grab a random chain
  randiter <- sample(x = unique(climeff$Iteration), size = 1) # grab a random MCMC iteration
  randyear <- sample(x = unique(intercept$yearid), size = 1) # grab a random year
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
  climcovs$ppt1TmeanSpr1 <- climcovs$ppt1*climcovs$TmeanSpr1
  climcovs$ppt2TmeanSpr2 <- climcovs$ppt2*climcovs$TmeanSpr2
  climcovs <- (climcovs - clim_scalers$means) / clim_scalers$sds
  
  cover[t] <- growFunc(N = cover[t-1], int = inttmp$value, 
                       slope = slopetmp$value, clims = tmpclim$value,
                       climcovs = climcovs, tau = tmptau$value) 
  setTxtProgressBar(pb, t)
}



####
####  Plot Simulation Results
####
par(mfrow=c(1,2))
plot(c(1:tsims), cover*100, type="l",ylab="Cover (%)",xlab="Time",main=do_species)
boxplot(cover*100, outline=FALSE, ylab="Cover (%)",main=do_species)
abline(h = mean(subset(allD, Species==do_species)[,"percCover"]*100), col="red")
text(x=1,y=mean(subset(allD, Species==do_species)[,"percCover"]*100),
     labels = "observed mean", col="red")

##  Print out results vs. observations
median(cover*100)
mean(subset(allD, Species==do_species)[,"percCover"]*100)
