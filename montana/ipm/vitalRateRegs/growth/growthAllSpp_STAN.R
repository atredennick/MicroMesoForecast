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
  vector[N] W; // crowding vector
}
parameters{
  real a_mu;
  real a[Yrs];
  real<lower=0> b1_mu;
  real<lower=0> b1[Yrs];
  vector[Covs] b2;
  real w;
  real gint[G];
  //real tau;
  //real tauSize;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sig_mod;
  real<lower=0> sig_G;
  real<lower=0> sig_res;
}
transformed parameters{
  real mu[N];
  //real tau2[N];
  //real tau3[N];
  vector[N] climEff;
  climEff <- C*b2;
  for(n in 1:N){
    mu[n] <- a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + w*W[n] + climEff[n];
    //tau2[n] <- tau*exp(tauSize*mu[n]); 
    //tau3[n] <- fmax(tau2[n],10000000);  
  }
}
model{
  // Priors
  a_mu ~ normal(0,1000);
  w ~ normal(0,1000);
  b1_mu ~ normal(0,1000);
  //tau ~ normal(0,1000);
  //tauSize ~ normal(0,1000);
  sig_a ~ uniform(0,1000);
  sig_b1 ~ uniform(0,1000);
  sig_mod ~ uniform(0,1000);
  sig_G ~ uniform(0,1000);
  for(g in 1:G)
      gint[g] ~ normal(0, sig_G);
  for(c in 1:Covs)
    b2[c] ~ uniform(-10,10);
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
    b1[y] ~ normal(b1_mu, sig_b1);
  }

  // Likelihood
  Y ~ normal(mu, sig_mod);
}
"

# rng_seed <- 123

## Loop through species and fit model
library(lme4)
big_list <- list()
for(do_species in sppList){
  growD <- subset(growD_all, species==do_species)
  outlm=lmer(log(area.t1)~log(area.t0)+W+pptLag+ppt1+TmeanSpr1+ 
             ppt2+TmeanSpr2+
             ppt1:TmeanSpr1+ppt2:TmeanSpr2+
             (log(area.t0)|year),data=subset(growD)) 
  summary(outlm)
  
  clim_covs <- growD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
  clim_covs$inter1 <- clim_covs$ppt1*clim_covs$TmeanSpr1
  clim_covs$inter2 <- clim_covs$ppt2*clim_covs$TmeanSpr2
  groups <- as.numeric(growD$Group)
  
  G <- length(unique(growD$Group))
  Yrs <- length(unique(growD$year))
  init_fun <- function(){list(a_mu=rnorm(1,0,1), a=rnorm(Yrs,0,1),
                              b1_mu=rnorm(1,0,1), b1=rnorm(Yrs,0,1),
                              b2=rnorm(length(clim_covs),0,1), gint=rnorm(G,0,1),
                              w=rnorm(1,0,1), tau=rnorm(1,0,1), tauSize=rnorm(1,0,1),
                              sig_a=runif(1,0,10), sig_b1=runif(1,0,10), sig_G=runif(1,0,10))}
  
  
  datalist <- list(N=nrow(growD), Yrs=Yrs, yid=(growD$year-31),
                   Covs=length(clim_covs), Y=log(growD$area.t1), X=log(growD$area.t0),
                   C=clim_covs, W=growD$W, G=G, gid=groups)
  pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
         "w", "gint", "sig_mod")
  mcmc_samples <- stan(model_code=model_string, data=datalist,
                       pars=pars,
                       chains=0)
#   traceplot(mcmc_samples)
  chain2 <- stan(fit=mcmc_samples, data=datalist, pars=pars,
                 chains=1, iter=2000, warmup=1000)
  plot(chain2)
  traceplot(chain2)
  print(chain2)
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







