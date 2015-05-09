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
                   area=NA,
                   survives=NA,
                   age=NA,
                   distEdgeMin=NA,
                   allEdge=NA,
                   species=NA)

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/edited/survD.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/survD.csv", sep=""))
    sppD$species <- doSpp 
  }
  outD <- rbind(outD, sppD)
}

survD <- outD[2:nrow(outD),]

climD <- read.csv("../../../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD[,clim_vars] <- scale(climD[,clim_vars], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900
survD <- merge(survD,climD)
survD$Group=as.factor(substr(survD$quad,1,1))

# Read in previously estimated crowding indices
c1 <- read.csv("BOGRsurvCrowding.csv")[,2:3]
c1$species <- sppList[1]
c2 <- read.csv("HECOsurvCrowding.csv")[,2:3]
c2$species <- sppList[2]
c3 <- read.csv("PASMsurvCrowding.csv")[,2:3]
c3$species <- sppList[3]
c4 <- read.csv("POSEsurvCrowding.csv")[,2:3]
c4$species <- sppList[4]
crowd <- rbind(c1,c2,c3,c4)

# Merge crowding and growth data
survD_all <- merge(survD, crowd, by=c("species", "X"))

#try glm
# fit final mixed effect model: based on what?
survD <- subset(survD_all, species=="BOGR")
# library(lme4)
# outlm=glmer(survives~log(area)+W+pptLag+ppt1+TmeanSpr1+ 
#            ppt2+TmeanSpr2+
#            ppt1:TmeanSpr1+ppt2:TmeanSpr2+
#            (1|Group) + (log(area)|year),data=subset(survD, species=="BOGR"),
#            family=binomial) 
# summary(outlm)

model_string <- "
data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> Covs; // climate covariates
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0,upper=1> Y[N]; // observation vector
  matrix[N,Covs] C; // climate matrix
  vector[N] X; // size vector
  vector[N] W; // crowding vector
}
parameters{
  real a_mu;
  vector[Yrs] a;
  real b1_mu;
  vector[Yrs] b1;
  vector[Covs] b2;
  real w;
  vector[G] gint;
  real<lower=0> sig_a;
  real<lower=0> sig_b1;
  real<lower=0> sig_G;
}
transformed parameters{
  real mu[N];
  vector[N] climEff;
  climEff <- C*b2;
  for(n in 1:N){
    mu[n] <- inv_logit(a[yid[n]] + gint[gid[n]] + b1[yid[n]]*X[n] + w*W[n] + climEff[n]);
  }
}
model{
  // Priors
  a_mu ~ uniform(-300,300);
  w ~ uniform(-100,100);
  b1_mu ~ uniform(-100,100);
  sig_a ~ cauchy(0,5);
  sig_b1 ~ cauchy(0,5);
  sig_G ~ cauchy(0,5);
  //for(g in 1:G)
      gint ~ normal(0, sig_G);
  //for(c in 1:Covs)
    b2 ~ uniform(-10,10);
  //for(y in 1:Yrs){
    a ~ normal(a_mu, sig_a);
    b1 ~ normal(b1_mu, sig_b1);
  //}

  // Likelihood
  //Y ~ bernoulli_logit(mu);
    Y ~ binomial(1,mu);
}
"

# rng_seed <- 123

## Loop through species and fit model
# library(lme4)
big_list <- list()

## Compile model outside of loop
survD <- subset(survD_all, species==sppList[1])
clim_covs <- survD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
clim_covs$inter1 <- clim_covs$ppt1*clim_covs$TmeanSpr1
clim_covs$inter2 <- clim_covs$ppt2*clim_covs$TmeanSpr2
groups <- as.numeric(survD$Group)
G <- length(unique(survD$Group))
Yrs <- length(unique(survD$year))

datalist <- list(N=nrow(survD), Yrs=Yrs, yid=(survD$year-31),
                 Covs=length(clim_covs), Y=survD$survives, X=log(survD$area),
                 C=clim_covs, W=survD$W, G=G, gid=groups)
pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
       "w", "gint")
mcmc_samples <- stan(model_code=model_string, data=datalist,
                     pars=pars, chains=0)


## Loop through and fit each species' model
for(do_species in sppList){
  survD <- subset(survD_all, species==do_species)
  
  clim_covs <- survD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
  clim_covs$inter1 <- clim_covs$ppt1*clim_covs$TmeanSpr1
  clim_covs$inter2 <- clim_covs$ppt2*clim_covs$TmeanSpr2
  groups <- as.numeric(survD$Group)
  G <- length(unique(survD$Group))
  Yrs <- length(unique(survD$year))
  
  ## Set reasonable initial values for three chains
  inits <- list()
  inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                     gint=rep(0,G), w=0, sig_b1=0.5, sig_a=0.5,
                     sig_G=0.5, b2=rep(0,length(clim_covs)))
  inits[[2]] <- list(a_mu=1, a=rep(1,Yrs), b1_mu=1, b1=rep(1,Yrs),
                     gint=rep(1,G), w=0.5, sig_b1=1, sig_a=1,
                     sig_G=1, b2=rep(1,length(clim_covs)))
  inits[[3]] <- list(a_mu=0.5, a=rep(0.5,Yrs), b1_mu=0.5, b1=rep(0.5,Yrs),
                     gint=rep(-1,G), w=-0.5, sig_b1=0.1, sig_a=0.1,
                     sig_G=0.1, b2=rep(-1,length(clim_covs)))
  
  datalist <- list(N=nrow(survD), Yrs=Yrs, yid=(survD$year-31),
                   Covs=length(clim_covs), Y=survD$survives, X=log(survD$area),
                   C=clim_covs, W=survD$W, G=G, gid=groups)
  pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
         "w", "gint")
  
  fitted <- stan(fit=mcmc_samples, data=datalist, pars=pars,
                 chains=3, iter=1000, warmup=250, init=inits)
  
  big_list[[do_species]] <- fitted
} # end species loop

saveRDS(big_list, "survival_stanfits_allspp.RDS")



