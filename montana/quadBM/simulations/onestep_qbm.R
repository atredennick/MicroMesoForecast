#Quad-Based Model simulations: one-step ahead forecasts

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

##  Set some global parameters
do_species <- "HECO"
outfile <- paste("./validation_results/onestep_", do_species, ".RDS", sep="")
NumberSimsPerYear <- 10

library(reshape2)
library(plyr)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(EnvStats)
library(rstan)
library(ggmcmc)

#bring in data
#bring in climate data
climD <- read.csv("../../weather/Climate.csv")
climD[2:6] <- scale(climD[2:6], center = TRUE, scale = TRUE)

allD <- read.csv("../../speciesData/quadAllCover.csv")
allD <- allD[,2:ncol(allD)] #get rid of X ID column
allD$percCover <- allD$totCover/10000
sppList <- as.character(unique(allD$Species))
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
allD <- backD[2:nrow(backD),]

##  Load vital rate parameters
fitlong <- readRDS(paste("../vitalRateRegressions/truncNormModel/popgrowth_stanmcmc_", 
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

# Group effects
group_off <- fitthin[grep("gint", fitthin$Parameter),]
group_off$groupid <- substr(group_off$Parameter, 6, length(group_off$Parameter))
group_off$groupid <- unlist(strsplit(group_off$groupid, split=']'))

# Lognormal sigma (called tau here)
tau <- fitthin[grep("tau", fitthin$Parameter),]

##  Define population growth function
growFunc <- function(N, int, slope, clims, climcovs, tau){
  mu <- mu <- int+slope*log(N)+sum(clims[1:7]*climcovs)+sum(clims[8:12]*log(N)*climcovs[1:5])
  newN <- rlnormTrunc(1, meanlog = mu, sdlog = tau, min = 0, max = 1)
  return(newN)
}


##  Run one-step forecasts
sppD <- subset(allD, Species==do_species)
quadnames <- as.data.frame(unique(allD$group))
colnames(quadnames) <- "group"
quadnames$groupnums <- as.numeric(as.factor(quadnames$group))
nSim <- NumberSimsPerYear
yearsN <- length(unique(sppD$year))
years <- as.numeric(unique(sppD$year)+1900)
yearsID <- unique(sppD$year)

output <- data.frame(year=NA, lagcov=NA, obscov=NA, predcov=NA, 
                     group=NA, rep=NA, species=NA)
for(do_year in 1:yearsN){
  yr_data <- subset(sppD, year==yearsID[do_year])
  quads <- as.numeric(as.factor(yr_data$group))
  quad_names <- yr_data$group
  climcovs <- yr_data[1,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
  climcovs$inter1 <- climcovs$ppt1*climcovs$TmeanSpr1
  climcovs$inter2 <- climcovs$ppt2*climcovs$TmeanSpr2
  for(i in 1:nrow(yr_data)){
    ntmp <- yr_data[i,"percLagCover"]
    gtmp <- yr_data[i,"group"]
    gnum <- as.numeric(subset(quadnames, group==gtmp)["groupnums"])
    for(do_sim in 1:nSim){
      randchain <- sample(x = climeff$Chain, size = 1)
      randiter <- sample(x = climeff$Iteration, size = 1)
      inttmp <- subset(intercept, Chain==randchain & 
                         Iteration==randiter &
                         yearid==do_year)
      grouptmp <- subset(group_off, Chain==randchain & 
                           Iteration==randiter &
                           groupid==gnum)
      slopetmp <- subset(coveff, Chain==randchain & 
                           Iteration==randiter &
                           yearid==do_year)
      tmpclim <- subset(climeff, Chain==randchain & 
                          Iteration==randiter)
      tmptau <- subset(tau, Chain==randchain & 
                         Iteration==randiter)
      
      # Get predicted cover
      covertmp <- growFunc(N = ntmp, int = (inttmp$value+grouptmp$value), 
                           slope = slopetmp$value, clims = tmpclim$value,
                           climcovs = climcovs, tau = tmptau$value)
      
      # Store in output dataframe
      tmpout <- data.frame(year=yearsID[do_year], lagcov=ntmp, 
                           obscov=yr_data[i,"percCover"], predcov=covertmp,
                           group=gtmp, rep=do_sim, species=do_species)
      output <- rbind(output, tmpout)
      print(paste("Done with rep", do_sim, "of observation", i, "for year", do_year))
    }#end quad-year rep loop
  }#end within-year observation loop
}#end year loop

output <- output[2:nrow(output),]
saveRDS(output, outfile)
# plot(output$obscov*100, output$predcov*100)
# abline(0,1, col="red")
