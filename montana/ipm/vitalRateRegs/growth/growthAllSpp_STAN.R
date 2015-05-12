#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rstan)
library(parallel)

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
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/edited/growDnoNA.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/growDnoNA.csv", sep=""))
    sppD$species <- doSpp 
  }
  outD <- rbind(outD, sppD)
}

growD <- outD[2:nrow(outD),]

climD <- read.csv("../../../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD[,clim_vars] <- scale(climD[,clim_vars], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900
growD <- merge(growD,climD)
growD$Group=as.factor(substr(growD$quad,1,1))

# Read in previously estimated crowding indices
c1 <- read.csv("BOGRgrowthCrowding.csv")[,2:3]
c1$species <- sppList[1]
c2 <- read.csv("HECOgrowthCrowding.csv")[,2:3]
c2$species <- sppList[2]
c3 <- read.csv("PASMgrowthCrowding.csv")[,2:3]
c3$species <- sppList[3]
c4 <- read.csv("POSEgrowthCrowding.csv")[,2:3]
c4$species <- sppList[4]
crowd <- rbind(c1,c2,c3,c4)

# crowd <- rbind(read.csv("BOGRgrowthCrowding.csv")[,2:3], 
#            read.csv("HECOgrowthCrowding.csv")[,2:3],
#            read.csv("PASMgrowthCrowding.csv")[,2:3],
#            read.csv("POSEgrowthCrowding.csv")[,2:3])
# crowd <- read.csv("HECOgrowthCrowding.csv")[,2:3] #ignore first column of rownames

# Merge crowding and growth data
growD_all <- merge(growD, crowd, by=c("species", "X"))

#try glm
# fit final mixed effect model: based on what?
# growD <- subset(growD_all, species=="HECO")
# library(lme4)
# outlm=lmer(log(area.t1)~log(area.t0)+W+pptLag+ppt1+TmeanSpr1+ 
#            ppt2+TmeanSpr2+
#            ppt1:TmeanSpr1+ppt2:TmeanSpr2+
#            (log(area.t0)|year),data=subset(growD, species=="HECO")) 
# summary(outlm)

model_string <- "
data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> Covs; // climate covariates
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  vector[N] Y; // observation vector
  matrix[N,Covs] C; // climate matrix
  vector[N] X; // size vector
  matrix[N,2] W; // crowding matrix
}
parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  vector[Covs] b2;
  vector[2] w;
  real gint[G];
  real tau;
  real tauSize;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sig_G;
}
transformed parameters{
  real mu[N];
  real<lower=0> tau2[N];
  vector[N] climEff;
  vector[N] crowdEff;
  climEff <- C*b2;
  crowdEff <- W*w;
  for(n in 1:N){
    mu[n] <- a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + crowdEff[n] + climEff[n];
    tau2[n] <- fmax(tau*exp(tauSize*mu[n]), 0.0000001);  
  }
}
model{
  // Priors
  a_mu ~ normal(0,1000);
  w ~ uniform(-100,100);
  b1_mu ~ normal(0,1000);
  tau ~ normal(0,1000);
  tauSize ~ normal(0,1000);
  sig_a ~ uniform(0,1000);
  sig_b1 ~ uniform(0,1000);
  sig_G ~ uniform(0,1000);
  for(g in 1:G)
      gint[g] ~ normal(0, sig_G);
  for(c in 1:Covs)
    b2[c] ~ normal(0,1000);
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
    b1[y] ~ normal(b1_mu, sig_b1);
  }

  // Likelihood
  Y ~ normal(mu, tau2);
}
"

# rng_seed <- 123

## Loop through species and fit model
# library(lme4)
big_list <- list()

## Compile model outside of loop
growD <- subset(growD_all, species==sppList[1])
clim_covs <- growD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
clim_covs$inter1 <- clim_covs$ppt1*clim_covs$TmeanSpr1
clim_covs$inter2 <- clim_covs$ppt2*clim_covs$TmeanSpr2
clim_covs$sizepptLag <- clim_covs$pptLag*log(growD$area.t0)
clim_covs$sizeppt1 <- clim_covs$ppt1*log(growD$area.t0)
clim_covs$sizeppt2 <- clim_covs$ppt2*log(growD$area.t0)
clim_covs$sizetemp1 <- clim_covs$TmeanSpr1*log(growD$area.t0)
clim_covs$sizetemp2 <- clim_covs$TmeanSpr2*log(growD$area.t0)
groups <- as.numeric(growD$Group)
G <- length(unique(growD$Group))
Yrs <- length(unique(growD$year))
W <- cbind(growD$W, growD$W*log(growD$area.t0))

datalist <- list(N=nrow(growD), Yrs=Yrs, yid=(growD$year-31),
                 Covs=length(clim_covs), Y=log(growD$area.t1), X=log(growD$area.t0),
                 C=clim_covs, W=W, G=G, gid=groups)
pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
       "w", "gint", "tau", "tauSize")
mcmc_samples <- stan(model_code=model_string, data=datalist,
                     pars=pars, chains=0)


## Loop through and fit each species' model
for(do_species in sppList){
  growD <- subset(growD_all, species==do_species)
  
  clim_covs <- growD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
  clim_covs$inter1 <- clim_covs$ppt1*clim_covs$TmeanSpr1
  clim_covs$inter2 <- clim_covs$ppt2*clim_covs$TmeanSpr2
  clim_covs$sizepptLag <- clim_covs$pptLag*log(growD$area.t0)
  clim_covs$sizeppt1 <- clim_covs$ppt1*log(growD$area.t0)
  clim_covs$sizeppt2 <- clim_covs$ppt2*log(growD$area.t0)
  clim_covs$sizetemp1 <- clim_covs$TmeanSpr1*log(growD$area.t0)
  clim_covs$sizetemp2 <- clim_covs$TmeanSpr2*log(growD$area.t0)
  groups <- as.numeric(growD$Group)
  G <- length(unique(growD$Group))
  Yrs <- length(unique(growD$year))
  W <- cbind(growD$W, growD$W*log(growD$area.t0))
  
  ## Set reasonable initial values for three chains
  inits <- list()
  inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                     gint=rep(0,G), w=c(0,0), sig_b1=0.5, sig_a=0.5, tau=0.5, tauSize=0.5,
                     sig_G=0.5, b2=rep(0,length(clim_covs)))
  inits[[2]] <- list(a_mu=1, a=rep(1,Yrs), b1_mu=1, b1=rep(1,Yrs),
                     gint=rep(1,G), w=c(0.5,0.5), sig_b1=1, sig_a=1, tau=1, tauSize=1,
                     sig_G=1, b2=rep(1,length(clim_covs)))
  inits[[3]] <- list(a_mu=0.5, a=rep(0.5,Yrs), b1_mu=0.5, b1=rep(0.5,Yrs),
                     gint=rep(-1,G), w=c(-0.5,-0.5), sig_b1=0.1, sig_a=0.1, tau=0.1, tauSize=0.1,
                     sig_G=0.1, b2=rep(-1,length(clim_covs)))
  
  datalist <- list(N=nrow(growD), Yrs=Yrs, yid=(growD$year-31),
                   Covs=length(clim_covs), Y=log(growD$area.t1), X=log(growD$area.t0),
                   C=clim_covs, W=W, G=G, gid=groups)
  pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
         "w", "gint", "tau", "tauSize")
  
  fitted <- stan(fit=mcmc_samples, data=datalist, pars=pars,
                 chains=3, iter=2000, warmup=1000, init=inits)
  
  big_list[[do_species]] <- fitted
} # end species loop

saveRDS(big_list, "growth_stanfits_allspp.RDS")



