#Quad-Based Model simulations

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))
yearvec <- readRDS("../../ipm/simulations/random_year_effects_sequence.rds")
climvec <- readRDS("../../ipm/simulations/climate_year_sequence.rds")

# do_species <- "BOGR"
tsims <- 2500
burn.in <- 500
perc_change <- 0.01

library(reshape2)
library(plyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(EnvStats)

#bring in data
allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
head(scale(allD$percCover, center=TRUE, scale=TRUE))
sppList <- as.character(unique(allD$Species))

#perturb climate data
climD <- read.csv("../../weather/Climate.csv")
climScale <- scale(climD[2:6], center = TRUE, scale = TRUE)
climAvg <- apply(X = climD, MARGIN = 2, FUN = mean)
climSD <- apply(X = climD, MARGIN = 2, FUN = sd)

pptVars=grep("ppt",names(climD))
tmp1=perc_change*colMeans(climD)
tmp1=matrix(tmp1,NROW(climD),length(tmp1),byrow=T)
climD[,pptVars]=climD[,pptVars]+tmp1[,pptVars]

climD[2] <- (climD[2] - climAvg[2])/climSD[2]
climD[3] <- (climD[3] - climAvg[3])/climSD[3]
climD[4] <- (climD[4] - climAvg[4])/climSD[4]
climD[5] <- (climD[5] - climAvg[5])/climSD[5]
climD[6] <- (climD[6] - climAvg[6])/climSD[6]


##  Loop through species
for(do_species in sppList){
  outfile <- paste("./results/climatechange_meanparams", do_species, "_qbm_cover_pptChange.RDS", sep="")
  
  ##  Load vital rate parameters
  fitlong <- readRDS(paste("../vitalRateRegressions/truncNormModel/popgrowth_stanmcmc_", 
                           do_species, ".RDS", sep=""))
  fitthin <- ddply(fitlong, .(Parameter), summarise,
                   value = mean(value))
  
  ##  Break up MCMC into regression components
  # Climate effects
  climeff <- fitthin[grep("b2", fitthin$Parameter),]
  climeff$id <- substr(climeff$Parameter, 4, length(climeff$Parameter))
  climeff$id <- unlist(strsplit(climeff$id, split=']'))
  climeff <- climeff[with(climeff, order(as.numeric(id))),]
  
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
  
  ##  Run simulations
  outD <- data.frame(cover=NA, species=NA, year=NA)
  cover <- numeric(tsims)
  cover[1] <- 0.01
  pb <- txtProgressBar(min=2, max=tsims, char="+", style=3, width=65)
  for(t in 2:tsims){
    inttmp <- intercept
    slopetmp <- coveff
    tmpclim <- climeff
    tmptau <- tau
    climyear <- climvec[t] - min(climvec)+1
    climcovs <- climD[climyear,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
    climcovs$inter1 <- climcovs$ppt1*climcovs$TmeanSpr1
    climcovs$inter2 <- climcovs$ppt2*climcovs$TmeanSpr2
    cover[t] <- growFunc(N = cover[t-1], int = inttmp$value, 
                         slope = slopetmp$value, clims = tmpclim$value,
                         climcovs = climcovs, tau = tmptau$value) 
    setTxtProgressBar(pb, t)
  }#end simulation loop
  # Save the output
  covd <- as.data.frame(cover[burn.in+1:tsims])
  covd$species <- do_species
  covd$climsim <- "ppt"
  colnames(covd)[1] <- "cover"
  saveRDS(covd, outfile)
}#end species loop

