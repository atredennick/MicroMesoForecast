##  Script to cycle through leave-one-year-out fits
##    of growth regressions via MCMC. These fits are then
##    used in IPMs to forecast the year missing from the
##    fitting procedure.

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  4-22-2015

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

## Set do_year for validation from command line prompt
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
leave_out_year <- as.numeric(myargument)
# leave_out_year <- 33

####
####  Load libraries ----------------------------------------------------------
####
library(rjags)
library(coda) 
load.module("dic")


####
#### Read in data by species and make one long data frame ---------------------
####
sppList=sort(c("BOGR","HECO","PASM","POSE"))    # Set list of species codes
outD <- data.frame(X=NA,
                   quad=NA,
                   year=NA,
                   trackID=NA,
                   area.t1=NA,
                   area.t0=NA,
                   age=NA,
                   species=NA)

# Begin looping through species list to load data
for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  if(doSpp == "BOGR"){ # IF BOGR, load the edited data
    sppD <- read.csv(paste("../../../speciesData/", 
                           doSpp, "/edited/growDnoNA.csv", sep=""))
  } else{ # Otherwise use normal path
    sppD <- read.csv(paste("../../../speciesData/", 
                           doSpp, "/growDnoNA.csv", sep=""))
  }
  # Grab just the colomns we want
  tmp_colnames <- which(colnames(outD)!="species")
  cols_to_keep <- which(colnames(sppD) %in% colnames(outD)[tmp_colnames])
  sppD <- sppD[,cols_to_keep]
  sppD$species <- doSpp
  outD <- rbind(outD, sppD)
}

growD <- outD[2:nrow(outD),]

climD <- read.csv("../../../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD[,clim_vars] <- scale(climD[,clim_vars], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900
growD <- merge(growD,climD)
growD$Group=as.factor(substr(growD$quad,1,1))

##  TODO: implement crowding via annuli indexing
# Read in previously estimated crowding indices
# crowd <- c(read.csv("BOGRgrowthCrowding.csv")[,2], 
#            read.csv("HECOgrowthCrowding.csv")[,2],
#            read.csv("PASMgrowthCrowding.csv")[,2],
#            read.csv("POSEgrowthCrowding.csv")[,2])
##  For checking script, include fake crowding indices = 0
crowd <- rep(x = 0, times = nrow(growD))

####
#### Set up data for JAGS model -----------------------------------------------
####
# Take out a year (set programmatically)
yr_data <- subset(growD, year!=leave_out_year)

dataJ <- list(Y = log(yr_data$area.t1),
              X = log(yr_data$area.t0),
              nObs = nrow(yr_data),
              grp = as.numeric(yr_data$Group),
              nGrp = length(unique(yr_data$Group)),
              spp = as.numeric(as.factor(yr_data$species)),
              nSpp = length(unique(yr_data$species)),
              yrs = as.numeric(as.factor((yr_data$year)-31)),
              nYrs = length(unique(yr_data$year)),
              crowd = crowd,
              TmeanSpr1 = yr_data$TmeanSpr1,
              TmeanSpr2 = yr_data$TmeanSpr2,
              ppt1 = yr_data$ppt1,
              ppt2 = yr_data$ppt2,
              pptlag = yr_data$pptLag)

# Set up some initial values for the MCMC chains (3 chains)
nyrs <- length(unique(yr_data$year))
nspp <- length(sppList)
ngrp <- length(unique(yr_data$Group))
inits=list(1)
inits[[1]]=list(intercept=rep(1,nspp), intYr=matrix(1, ncol=nyrs, nrow=nspp),
                intG=matrix(0.1,ncol=ngrp, nrow=nspp), betaSpp=rep(1,nspp), betaVar=rep(1,nspp),
                beta=matrix(1,ncol=nyrs,nrow=nspp), nb=rep(-0.2,nspp), intVarG=rep(2,nspp),
                tau=rep(1,nspp), tauSize=rep(1,nspp), intVaryY=rep(1,nspp),
                temp1Mu=1, temp2Mu=1, rain1Mu=1, rain2Mu=1, rainlagMu=1,
                temp1Var=1, temp2Var=1, rain1Var=1, rain2Var=1, rainlagVar=1,
                temp1=rep(1,nspp), temp2=rep(1,nspp), rain1=rep(1,nspp), rain2=rep(1,nspp), rainlag=rep(1,nspp))

inits[[2]]=list(intercept=rep(-1,nspp), intYr=matrix(-1, ncol=nyrs, nrow=nspp),
                intG=matrix(-0.1,ncol=ngrp, nrow=nspp), betaSpp=rep(-1,nspp), betaVar=rep(0.5,nspp),
                beta=matrix(-1,ncol=nyrs,nrow=nspp), nb=rep(0.2,nspp), intVarG=rep(0.5,nspp),
                tau=rep(0.5,nspp), tauSize=rep(0.5,nspp), intVaryY=rep(0.5,nspp),
                temp1Mu=2, temp2Mu=2, rain1Mu=2, rain2Mu=2, rainlagMu=2,
                temp1Var=2, temp2Var=2, rain1Var=2, rain2Var=2, rainlagVar=1,
                temp1=rep(2,nspp), temp2=rep(2,nspp), rain1=rep(2,nspp), rain2=rep(2,nspp), rainlag=rep(2,nspp))

inits[[3]]=list(intercept=rep(0.1,nspp), intYr=matrix(0.1, ncol=nyrs, nrow=nspp),
                intG=matrix(0.5,ncol=ngrp, nrow=nspp), betaSpp=rep(0.1,nspp), betaVar=rep(0.1,nspp),
                beta=matrix(0.1,ncol=nyrs,nrow=nspp), nb=rep(-0.5,nspp), intVarG=rep(0.2,nspp),
                tau=rep(0.1,nspp), tauSize=rep(0.1,nspp), intVaryY=rep(0.1,nspp),
                temp1Mu=0.1, temp2Mu=0.1, rain1Mu=0.1, rain2Mu=0.1, rainlagMu=0.1,
                temp1Var=0.1, temp2Var=0.1, rain1Var=0.1, rain2Var=0.1, rainlagVar=0.1,
                temp1=rep(0.1,nspp), temp2=rep(0.1,nspp), rain1=rep(0.1,nspp), rain2=rep(0.1,nspp), rainlag=rep(0.1,nspp))


####
#### Run MCMC from JAGS ------------------------
####
iterations <- 50000
adapt <- 5000
mod <- jags.model("growthAllSpp_JAGS.R", data=dataJ, n.chains=length(inits), 
                  n.adapt=adapt, inits=inits)
update(mod, n.iter = (iterations*0.25))
out <- coda.samples(mod, c("intYr", "beta", "intG", "nb", "temp1", "temp2", 
                           "rain1", "rain2", "rainlag", "intercept", "betaSpp"),
                    n.iter=iterations, n.thin=10)

####
#### Convert to dataframe for export and get other summaries ------------------
####
gelmDiag <- gelman.diag(out)

outC <- rbind(out[[1]][(iterations-999):iterations,], 
              out[[2]][(iterations-999):iterations,], 
              out[[3]][(iterations-999):iterations,])

outStat <- as.data.frame(summary(out)$stat)
outQuant <- as.data.frame(summary(out)$quantile)

sppNames <- c(rep(sppList, 12+6+12+4+4))
outStat$species <- sppNames
outQuant$species <- sppNames

uniq_years <- unique(yr_data$year)
year_names <- c(rep(uniq_years, each=4), rep(NA,7*4),
                rep(uniq_years, each=4), rep(NA,7*4))
outStat$year <- year_names
outQuant$year <- year_names

saveRDS(outC, file = paste("growthParamsMCMC_", leave_out_year, ".rds", sep=""))
write.csv(gelmDiag[[1]], file=paste("growthGelman_", leave_out_year, ".csv", sep=""))
write.csv(outStat, file=paste("growthStats_", leave_out_year, ".csv", sep=""))
write.csv(outQuant, file=paste("growthQuants_", leave_out_year, ".csv", sep=""))


