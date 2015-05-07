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
  real Y[N]; // observation vector
  real C[N,Covs]; // climate matrix
  real X[N]; // size vector
  real W[N]; // crowding vector
}
parameters{
  real a_mu;
  real a[Yrs];
  real b1_mu;
  real b1[Yrs];
  real b2[Covs];
  real w;
  real gint[G];
  real tau;
  real tauSize;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  //real<lower=0> sig_mod;
  real<lower=0> sig_G;
}
transformed parameters{
  real mu[N];
  real tau2[N];
  real tau3[N];
  for(n in 1:N){
    mu[n] <- a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + w*W[n] + b2[1]*C[n,1] + b2[2]*C[n,2] + b2[3]*C[n,3] + b2[4]*C[n,4] + b2[5]*C[n,5] + b2[6]*C[n,6];
    tau2[n] <- 1/(tau*exp(tauSize*mu[n])); 
    tau3[n] <- fmax(tau2[n],0.00000001);  
  }
}
model{
  // Priors
  a_mu ~ normal(0,10);
  w ~ uniform(-5,5);
  b1_mu ~ normal(0,1);
  tau ~ normal(0,1000);
  tauSize ~ normal(0,100);
  sig_a ~ uniform(0,1000);
  sig_b1 ~ uniform(0,1000);
  //sig_mod ~ uniform(0,1000);
  sig_G ~ uniform(0,1000);
  for(g in 1:G)
      gint[g] ~ normal(0, sig_G);
  for(c in 1:Covs)
    b2[c] ~ normal(0,10);
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
    b1[y] ~ normal(b1_mu, sig_b1);
  }

  // Likelihood
  Y ~ normal(mu, tau3);
}
"

# rng_seed <- 123

## Loop through species and fit model
big_list <- list()
for(do_species in sppList){
  growD <- subset(growD_all, species==do_species)
  clim_covs <- growD[,c("ppt1","ppt2","TmeanSpr1","TmeanSpr2")]
  clim_covs$inter1 <- clim_covs$ppt1*clim_covs$TmeanSpr1
  clim_covs$inter2 <- clim_covs$ppt2*clim_covs$TmeanSpr2
  groups <- as.numeric(growD$Group)
  datalist <- list(N=nrow(growD), Yrs=length(unique(growD$year)), yid=(growD$year-31),
                   Covs=length(clim_covs), Y=log(growD$area.t1), X=log(growD$area.t0),
                   C=clim_covs, W=growD$W, G=length(unique(growD$Group)), gid=groups)
  
  mcmc_samples <- stan(model_code=model_string, data=datalist,
                       pars=c("a_mu", "a", "b1_mu",  "b1", "b2", 
                              "w", "gint", "tau", "tauSize"),
                       chains=3, iter=1500, warmup=500)
  big_list[[do_species]] <- mcmc_samples
} # end species loop

# clim_covs <- growD[,c("ppt1","ppt2","TmeanSpr1","TmeanSpr2")]
# clim_covs$inter1 <- clim_covs$ppt1*clim_covs$TmeanSpr1
# clim_covs$inter2 <- clim_covs$ppt2*clim_covs$TmeanSpr2
# groups <- as.numeric(growD$Group)
# datalist <- list(N=nrow(growD), Yrs=length(unique(growD$year)), yid=(growD$year-31),
#                  Covs=length(clim_covs), Y=log(growD$area.t1), X=log(growD$area.t0),
#                  C=clim_covs, W=growD$W, G=length(unique(growD$Group)), gid=groups)
# 
# mcmc_samples <- stan(model_code=model_string, data=datalist,
#                      pars=c("b1", "b2", "w", "gint", "tau", "tauSize"),
#                      chains=3, iter=1000, warmup=150)
# 
# traceplot(mcmc_samples)
# plot(mcmc_samples)
# print(mcmc_samples)

# library(ggmcmc)
# ggmcs <- ggs(mcmc_samples)







