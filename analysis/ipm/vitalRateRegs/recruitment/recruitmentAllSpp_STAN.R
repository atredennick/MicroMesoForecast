#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rstan)
library(parallel)
library(reshape2)

##  Read in lppd scores and selected prior climate stddevs
priors_df <- read.csv("../../../all_maxlppds.csv")
priors <- subset(priors_df, vital=="recruitment")

sppList=sort(c("BOGR","HECO","PASM","POSE"))

####
####  Load Climate Covariates and Recruitment Data
####
climD <- read.csv("../../../weather/Climate.csv") #on local machine
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD$year <- climD$year-1900



####
#### Read in data by species and make one long data frame -------------
####
for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/edited/recArea.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("../../../speciesData/", doSpp, "/recArea.csv", sep=""))
    sppD$species <- doSpp 
  }
  
  sppD$Group <- as.factor(substr(sppD$quad,1,1)) #add by Chengjin
  sppD <- sppD[,c("quad","year","NRquad","totParea","Group")]
  names(sppD)[3] <- paste("R.",sppList[spp],sep="")
  names(sppD)[4] <- paste("cov.",sppList[spp],sep="")
  if(spp==1){
    D <- sppD
  }else{
    D <- merge(D,sppD,all=T)
  }
}
D[is.na(D)]=0  # replace missing values 

# calculate mean cover by group and year
tmpD <- D[,c("quad","year","Group",paste("cov.",sppList,sep=""))]
tmpD <- aggregate(tmpD[,4:NCOL(tmpD)],
                  by=list("year"=tmpD$year,"Group"=tmpD$Group),FUN=mean)
names(tmpD)[3:NCOL(tmpD)] <- paste("Gcov.",sppList,sep="")
D <- merge(D,tmpD,all.x=T)




###if using square transform, there is no need of the following code
####################################################################
###here you need to check: both parents1 and parents2 are equal to 0 at the same time
parents1 <- as.matrix(D[,c(paste("cov.",sppList,sep=""))])/100 ##convert from absolute cover to [1,100] range
parents2 <- as.matrix(D[,c(paste("Gcov.",sppList,sep=""))])/100

##for species 1
tmp1L <- which(parents1[,1]==0) ##lcoal
tmp1G <- which(parents2[,1]==0) ##Group
tmp1 <- intersect(tmp1L,tmp1G)
##for species 2
tmp2L <- which(parents1[,2]==0)
tmp2G <- which(parents2[,2]==0)
tmp2 <- intersect(tmp2L,tmp2G)
##for species 3
tmp3L <- which(parents1[,3]==0)
tmp3G <- which(parents2[,3]==0)
tmp3 <- intersect(tmp3L,tmp3G)
##for species 4
tmp4L <- which(parents1[,4]==0)
tmp4G <- which(parents2[,4]==0)
tmp4 <- intersect(tmp4L,tmp4G)

tmp <- unique(c(tmp1,tmp2,tmp3,tmp4))

if(length(tmp)>0){
  parents1 <- parents1[-tmp,] ##remove them
  parents2 <-parents2[-tmp,] ##remove them
  y <- as.matrix(D[,c(paste("R.",sppList,sep=""))])[-tmp,] ##remove them  
  year <- D$year[-tmp] ##remove them
  Nyrs <- length(unique(D$year))
  N <- dim(D)[1]-length(tmp) ##reduce
  Nspp <- length(sppList)
  Group <- as.numeric(as.factor(D$Group))[-tmp] ##remove them ##first turn it as FACTOR, then to NUMERIC
  Ngroups <- length(unique(Group))
} else {
  y <- as.matrix(D[,c(paste("R.",sppList,sep=""))])
  year <- D$year
  Nyrs <- length(unique(D$year))
  N <- dim(D)[1]
  Nspp <- length(sppList)
  Group <- as.numeric(as.factor(D$Group)) ##first turn it as FACTOR, then to NUMERIC
  Ngroups <- length(unique(Group))
}

tmpY <- melt(y)
tmpP1 <- melt(parents1)
tmpP2 <- melt(parents2)
allD <- data.frame(species=tmpY$Var2,
                   year=rep(year,4),
                   group=rep(Group,4),
                   recruits=tmpY$value,
                   parents1=tmpP1$value,
                   parents2=tmpP2$value)

##  Add in climate data
allD <- merge(allD,climD)


# Stan model
model_string <- "
data{
  int<lower=0> N; // observations
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> Covs; // climate covariates
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> Y[N]; // observation vector
  matrix[N,Covs] C; // climate matrix
  vector[N] parents1; // crowding vector
  vector[N] parents2; // crowding vector
  real tauclim; // prior climate std dev
}
parameters{
  real a_mu;
  vector[Yrs] a;
  vector[Covs] b2;
  real dd;
  real gint[G];
  real<lower=0> sig_a;
  real<lower=0> theta;
  real<lower=0> sig_G;
  real<lower=0, upper=1> u;
}
transformed parameters{
  real mu[N];
  vector[N] climEff;
  vector[N] trueP1;
  vector[N] trueP2;
  vector[N] lambda;
  vector[N] q;
  climEff <- C*b2;
  for(n in 1:N){
    trueP1[n] <- parents1[n]*u + parents2[n]*(1-u);
    trueP2[n] <- sqrt(trueP1[n]);
    mu[n] <- exp(a[yid[n]] + gint[gid[n]] + dd*trueP2[n] + climEff[n]);
    lambda[n] <- trueP1[n]*mu[n];
    q[n] <- lambda[n]*theta;
  }
}
model{
  u ~ uniform(0,1);
  theta ~ uniform(0,100);
  a_mu ~ normal(0,1000);
  dd ~ uniform(-100,100);
  sig_a ~ cauchy(0,5);
  sig_G ~ cauchy(0,5);
  for(g in 1:G)
  gint[g] ~ normal(0, sig_G);
  for(y in 1:Yrs){
    a[y] ~ normal(a_mu, sig_a);
  }
  for(j in 1:Covs)
    b2[j] ~ normal(0, tauclim);

  // Likelihood
  Y ~ neg_binomial_2(q, theta);
}
"


## Compile model outside of loop
recD <- subset(allD, species==paste("R.",sppList[1],sep=""))
##  Create and scale interaction covariates
recD$ppt1TmeanSpr1 <- recD$ppt1*recD$TmeanSpr1
recD$ppt2TmeanSpr2 <- recD$ppt2*recD$TmeanSpr2
clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2")
clim_covs <- recD[,clim_vars_all]
clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
groups <- as.numeric(recD$group)
G <- length(unique(recD$group))
Yrs <- length(unique(recD$year))
yid <- as.numeric(as.factor(recD$year))

datalist <- list(N=nrow(recD), Yrs=Yrs, yid=yid,
                 Covs=ncol(clim_covs), Y=recD$recruits, C=clim_covs, 
                 parents1=recD$parents1, parents2=recD$parents2,
                 G=G, gid=groups, tauclim=0.1)
pars=c("a_mu", "a", "u", "theta", "b2",
       "dd", "gint")
mcmc_samples <- stan(model_code=model_string, data=datalist,
                     pars=pars, chains=0)


big_list <- list()
for(do_species in 1:length(sppList)){
  recD <- subset(allD, species==paste("R.",sppList[do_species],sep=""))
  ##  Create and scale interaction covariates
  recD$ppt1TmeanSpr1 <- recD$ppt1*recD$TmeanSpr1
  recD$ppt2TmeanSpr2 <- recD$ppt2*recD$TmeanSpr2
  clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2")
  clim_covs <- recD[,clim_vars_all]
  clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE)
  groups <- as.numeric(recD$group)
  G <- length(unique(recD$group))
  Yrs <- length(unique(recD$year))
  yid <- as.numeric(as.factor(recD$year))

  # Get climate prior std dev
  prior_stddev <- as.numeric(subset(priors, species==sppList[do_species])["prior_stdev"])
  
  datalist <- list(N=nrow(recD), Yrs=Yrs, yid=yid,
                   Covs=ncol(clim_covs), Y=recD$recruits, C=clim_covs, 
                   parents1=recD$parents1, parents2=recD$parents2,
                   G=G, gid=groups, tauclim=prior_stddev)
  pars=c("a_mu", "a", "u", "theta", "b2",
         "dd", "gint")
  
  inits=list()
  inits[[1]]=list(a=rep(4,Yrs), a_mu=1, sig_a=1,
                  gint=rep(0,G), sig_G=1, u=0.4,
                  dd=-1,theta=1, b2=rep(0,ncol(clim_covs))) 
  inits[[2]]=list(a=rep(0.5,Yrs), a_mu=0.2, sig_a=10,
                  gint=rep(0,G), sig_G=0.1,  u=0.7,
                  dd=-0.05,theta=1.5, b2=rep(0.5,ncol(clim_covs))) 
  inits[[3]]=list(a=rep(1,Yrs), a_mu=0.5, sig_a=5,
                  gint=rep(-0.1,G), sig_G=0.5,  u=0.5,
                  dd=-1,theta=1, b2=rep(-0.5,ncol(clim_covs))) 
  
#   fitted <- stan(fit=mcmc_samples, data=datalist, pars=pars,
#                  chains=3, iter=1000, warmup=150, init=inits)
  rng_seed <- 123
  sflist <-
    mclapply(1:3, mc.cores=3,
             function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                              seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                              iter=2000, warmup=1000, init=list(inits[[i]])))
  fit <- sflist2stanfit(sflist)
  big_list[[sppList[do_species]]] <- fit
}

saveRDS(big_list, "recruitment_stanfits.RDS")
