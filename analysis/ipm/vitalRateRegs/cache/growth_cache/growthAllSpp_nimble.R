#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library('nimble')
library('tidyr')

sppList=sort(c("BOGR","HECO","PASM","POSE"))

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

##then we moved some specific points:
tmp2<-which(growD$quad=="A12" & growD$year==44)
tmp3<-which(growD$quad=="B1"  & growD$year==44)
tmp41<-which(growD$quad=="E4" & growD$year==33) 
tmp42<-which(growD$quad=="E4" & growD$year==34) 
tmp43<-which(growD$quad=="E4" & growD$year==43)
tmp44<-which(growD$quad=="E4" & growD$year==44)
tmpONE<-c(tmp2,tmp3,tmp41,tmp42,tmp43,tmp44)
if(length(tmpONE)>0) growD<-growD[-tmpONE,]

climD <- read.csv("../../../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD[,clim_vars] <- scale(climD[,clim_vars], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900
growD <- merge(growD,climD)
growD$Group=as.factor(substr(growD$quad,1,1))

# Read in previously estimated crowding indices
crowd <- c(read.csv("BOGRgrowthCrowding.csv")[,2], 
           read.csv("HECOgrowthCrowding.csv")[,2],
           read.csv("PASMgrowthCrowding.csv")[,2],
           read.csv("POSEgrowthCrowding.csv")[,2])


####
#### Set up data for JAGS model ---------------------------
####
growD <- subset(growD, year < 34)
constants <- list(X = log(growD$area.t0),
              nObs = nrow(growD),
              grp = as.numeric(growD$Group),
              nGrp = length(unique(growD$Group)),
              spp = as.numeric(as.factor(growD$species)),
              nSpp = length(unique(growD$species)),
              yrs = (growD$year - 31),
              nYrs = length(unique(growD$year)),
              crowd = crowd,
              TmeanSpr1 = growD$TmeanSpr1,
              TmeanSpr2 = growD$TmeanSpr2,
              ppt1 = growD$ppt1,
              ppt2 = growD$ppt2,
              pptlag = growD$pptLag)
data <- list(Y = log(growD$area.t1))

# Set up some initial values for the MCMC chains (3 chains)
nyrs <- length(unique(growD$year))
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
####  BUGS model -----------------------------------------
####
growth_model <- nimbleCode({
  #process model and likelihood
  for(i in 1:nObs){
    mu[i] <- intYr[spp[i],yrs[i]] + intG[spp[i],grp[i]] + beta[spp[i],yrs[i]]*X[i] + 
      nb[spp[i]]*crowd[i] + temp1[spp[i]]*TmeanSpr1[i] + 
      temp2[spp[i]]*TmeanSpr2[i] + rain1[spp[i]]*ppt1[i] + 
      rain2[spp[i]]*ppt2[i] + rainlag[spp[i]]*pptlag[i]
    tau2[i] <- 1/(tau[spp[i]]*exp(tauSize[spp[i]]*mu[i])) 
    tau3[i] <- max(tau2[i],0.00000001)  
    Y[i] ~ dnorm(mu[i], tau3[spp[i]])
  }
  
  #priors
  for(s in 1:nSpp){
    tau[s] ~ dnorm(0,0.001)
    tauSize[s] ~ dnorm(0,0.001)
    betaSpp[s] ~ dnorm(0, 1e-6)
    nb[s] ~ dnorm(0, 1e-6)
    temp1[s] ~ dnorm(temp1Mu, temp1Var)
    temp2[s] ~ dnorm(temp2Mu, temp2Var)
    rain1[s] ~ dnorm(rain1Mu, rain1Var)
    rain2[s] ~ dnorm(rain2Mu, rain2Var)
    rainlag[s] ~ dnorm(rainlagMu, rainlagVar)
    intercept[s] ~ dnorm(0, 1e-6)
    intVaryY[s] ~ dgamma(0.001, 0.001)
    betaVar[s] ~ dgamma(0.001, 0.001)
    intVarG[s] ~ dgamma(2, 0.5) 
    for(y in 1:nYrs){
      intYr[s,y] ~ dnorm(intercept[s], intVaryY[s])
      beta[s,y] ~ dnorm(betaSpp[s], betaVar[s])
    }
    for(g in 1:nGrp){
      intG[s,g] ~ dnorm(0, intVarG[s])
    }
  }
  
  temp1Mu ~ dnorm(0,1e-6)
  temp2Mu ~ dnorm(0,1e-6)
  rain1Mu ~ dnorm(0,1e-6)
  rain2Mu ~ dnorm(0,1e-6)
  rainlagMu ~ dnorm(0,1e-6)
  temp1Var ~ dgamma(0.001, 0.001)
  temp2Var ~ dgamma(0.001, 0.001)
  rain1Var ~ dgamma(0.001, 0.001)
  rain2Var ~ dgamma(0.001, 0.001)
  rainlagVar ~ dgamma(0.001, 0.001)
})

####
####  Run NIMBLE MCMC ------------------------------------
####
climate_covariates <- c("temp1", "temp2", "rain1", "rain2", "rainlag")

ptm <- proc.time() #start timer

Rmodel <- nimbleModel(code = growth_model, 
                      constants = constants, 
                      data = data,
                      inits = inits[[1]])
mcmcspec <- configureMCMC(Rmodel, print=TRUE, thin=100)
mcmcspec$addSampler('RW_block', list(targetNodes = climate_covariates,
                                     adaptInterval = 100))
mcmcspec$addMonitors(c('beta', 'intG', climate_covariates))
Rmcmc <- buildMCMC(mcmcspec)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)
Cmcmc$run(2500)
Cmcmc$run(50000)

proc.time() - ptm #get processing time
samples <- as.data.frame(as.matrix(Cmcmc$mvSamples))
dim(samples)
# head(samples)
long <- gather(samples)
apply(samples, 2, mean)
ggplot(long) + 
  geom_line(aes(seq_along(value), value)) + 
  facet_wrap(~key, scale='free')