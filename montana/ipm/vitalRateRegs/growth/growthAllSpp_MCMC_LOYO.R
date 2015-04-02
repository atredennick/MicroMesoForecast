#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rjags)
library(coda)
library(parallel) 
library(snowfall) 
library(rlecuyer) 
load.module("dic")

sppList=sort(c("BOGR","HECO","PASM","POSE"))
year_index_for_leave_out <- 1

## Set do_year for validation from command line prompt
# args <- commandArgs(trailingOnly = F)
# myargument <- args[length(args)]
# myargument <- sub("-","",myargument)
# do_year <- as.numeric(myargument)

####
#### Read in data by species and make one long data frame -------------
####
outD <- data.frame(X=NA,
                   quad=NA,
                   year=NA,
                   trackID=NA,
                   area.t1=NA,
                   area.t0=NA,
                   age=NA,
                   allEdge=NA,
                   distEdgeMin=NA,
                   species=NA)

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  sppD <- read.csv(paste("../../../speciesData/", doSpp, "/growDnoNA.csv", sep=""))
  sppD$species <- doSpp
  outD <- rbind(outD, sppD)
}

growD <- outD[2:nrow(outD),]

##  Leave one year out for fitting
all_years <- unique(growD$year)
year_to_leave_out <- all_years[year_index_for_leave_out]
growD <- subset(growD, year!=year_to_leave_out)

climD <- read.csv("../../../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD[,clim_vars] <- scale(climD[,clim_vars], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900
growD <- merge(growD,climD)
growD$Group=as.factor(substr(growD$quad,1,1))

source("get_crowding_loyo.R")


####
#### Set up data for JAGS model ---------------------------
####
#this is to get the sequencing right for indexing in the JAGS model
yrs <- as.numeric(as.factor((growD$year)-32))
dataJ <- list(Y = log(growD$area.t1),
              X = log(growD$area.t0),
              nObs = nrow(growD),
              grp = as.numeric(growD$Group),
              nGrp = length(unique(growD$Group)),
              spp = as.numeric(as.factor(growD$species)),
              nSpp = length(unique(growD$species)),
              yrs = yrs,
              nYrs = length(unique(growD$year)),
              crowd = crowd,
              TmeanSpr1 = growD$TmeanSpr1,
              TmeanSpr2 = growD$TmeanSpr2,
              ppt1 = growD$ppt1,
              ppt2 = growD$ppt2,
              pptlag = growD$pptLag)

# Set up some initial values for the MCMC chains (3 chains)
nyrs <- length(unique(yrs))
nspp <- length(sppList)
ngrp <- length(unique(growD$Group))
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
iterations <- 100
adapt <- 100
mod <- jags.model("growthAllSpp_JAGS.R", data=dataJ, n.chains=length(inits), 
                  n.adapt=adapt, inits=inits)
update(mod, n.iter = (iterations*0.25))
out <- coda.samples(mod, c("intYr", "beta", "intG", "nb", "temp1", "temp2", 
                           "rain1", "rain2", "rainlag", "intercept", "betaSpp"),
                    n.iter=iterations, n.thin=10)

####
#### Check for convergence ---------------------
####
gelmDiag <- gelman.diag(out)

####
#### Convert to dataframe for export and get other summaries ------------------
####
outC <- rbind(out[[1]][(iterations-999):iterations,], 
              out[[2]][(iterations-999):iterations,], 
              out[[3]][(iterations-999):iterations,])
saveRDS(outC, file = paste("./LOYO_results/growthParamsMCMC_", year_to_leave_out, ".rds", sep=""))
write.csv(gelmDiag[[1]], file=paste("./LOYO_results/growthGelman_", year_to_leave_out, ".csv"))


