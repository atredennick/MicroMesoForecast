#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
# Set working directory to location of this source file #
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rstan)
library(parallel)
library(reshape2)

## Set do_year for validation from command line prompt
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
do_species <- as.numeric(myargument)

sppList=sort(c("BOGR","HECO","PASM","POSE"))
do_species <- 2

####
#### Read in data by species and make one long data frame -------------
####
for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste("../../../../speciesData/", doSpp, "/edited/recArea.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste("../../../../speciesData/", doSpp, "/recArea.csv", sep=""))
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


climD <- read.csv("../../../../weather/Climate.csv")
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
climD[,clim_vars] <- scale(climD[,clim_vars], center = TRUE, scale = TRUE)
climD$year <- climD$year-1900
allD <- merge(allD,climD)


model_string <- "
data{
  int<lower=0> N; // observations
  int<lower=0> npreds;
  int<lower=0> Yrs; // years
  int<lower=0> yid[N]; // year id
  int<lower=0> Covs; // climate covariates
  int<lower=0> G; // groups
  int<lower=0> gid[N]; // group id
  int<lower=0> gid_out[npreds]; // group id
  int<lower=0> Y[N]; // observation vector
  int<lower=0> y_holdout[npreds]; // observation vector
  matrix[N,Covs] C; // climate matrix
  matrix[npreds,Covs] C_out; // climate matrix, holdout
  vector[N] parents1; // crowding vector
  vector[N] parents2; // crowding vector
  vector[npreds] parents1_out; // crowding vector
  vector[npreds] parents2_out; // crowding vector
  real<lower=0> tau; // prior variance
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
  //vector<lower=0>[Covs] lambda1;
  //real<lower=0> tau;
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
// Priors
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
  //lambda1 ~ cauchy(0, 1);
  //tau ~ cauchy(0, 1);
  for(j in 1:Covs)
    //b2[j] ~ normal(0, lambda1[j] * tau);
    b2[j] ~ normal(0, tau);

// Likelihood
Y ~ neg_binomial_2(q, theta);
}
generated quantities {
  real muhat[npreds];
  vector[npreds] climpred;
  vector[npreds] trueP1_pred;
  vector[npreds] trueP2_pred;
  vector[npreds] lambda_hat;
  vector[npreds] qpred;
  vector[npreds] lppd; // vector for computing log pointwise predictive density
  climpred <- C_out*b2;
  for(n in 1:npreds){
    trueP1_pred[n] <- parents1_out[n]*u + parents2_out[n]*(1-u);
    trueP2_pred[n] <- sqrt(trueP1_pred[n]);
    muhat[n] <- exp(a_mu + gint[gid_out[n]] + dd*trueP2_pred[n] + climpred[n]);
    lambda_hat[n] <- trueP1_pred[n]*muhat[n];
    qpred[n] <- lambda_hat[n]*theta;
    lppd[n] <- neg_binomial_2_log(y_holdout[n], qpred[n], theta);
  } 
}
"


## Compile model outside of loop
recD <- subset(allD, species==paste("R.",sppList[2],sep=""))
clim_covs <- recD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
groups <- as.numeric(recD$group)
G <- length(unique(recD$group))
Yrs <- length(unique(recD$year))

holdD <- subset(recD, year==32)
clim_holds <- holdD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
groups_out <- as.numeric(holdD$group)
Gout <- length(unique(holdD$group))
npreds <- nrow(holdD)
y_holdout <- holdD$recruits
parents1_out <- holdD$parents1
parents2_out <- holdD$parents2


datalist <- list(N=nrow(recD), Yrs=Yrs, yid=(recD$year-31),
                 Covs=length(clim_covs), Y=recD$recruits, C=clim_covs, 
                 parents1=recD$parents1, parents2=recD$parents2,
                 G=G, gid=groups, tau=(1000),
                 npreds=npreds, y_holdout=y_holdout, C_out=clim_holds,
                 parents1_out=parents1_out, parents2_out=parents2_out,
                 gid_out=groups_out, Gout=Gout)

pars=c("lppd")
mcmc_samples <- stan(model_code=model_string, data=datalist,
                     pars=pars, chains=0)

iterations <- 200
n.mcmc <- iterations/2
tau.set <- c(1000,500,100,50,10,5,1,0.1)
lambda.set <- 1/tau.set
# tau.set <- 1/lambda.set
nlam <- length(tau.set)
ncvs <- length(unique(recD$year))
year_vec <- unique(recD$year)

lppd <- rep(NA, length(year_vec))
lppd2 <- matrix(NA, nrow=length(year_vec), ncol=length(tau.set))

beta_cvs <- matrix(NA, nrow=ncvs, ncol=ncol(clim_covs))
beta_regs <- array(NA, c(ncvs, ncol(clim_covs), length(tau.set)))


recD_all <- recD
for(j in 1:length(tau.set)){
  for(s in 1:ncvs){
    recD <- subset(recD_all, year != year_vec[s])
    holdD <- subset(recD_all, year==year_vec[s])
    clim_covs <- recD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
    groups <- as.numeric(recD$group)
    G <- length(unique(recD$group))
    Yrs <- length(unique(recD$year))
    yid <- as.numeric(as.factor(recD$year))
    clim_holds <- holdD[,c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")]
    groups_out <- as.numeric(holdD$group)
    Gout <- length(unique(holdD$group))
    npreds <- nrow(holdD)
    y_holdout <- holdD$recruits
    parents1_out <- holdD$parents1
    parents2_out <- holdD$parents2
    
    
    datalist <- list(N=nrow(recD), Yrs=Yrs, yid=yid,
                     Covs=length(clim_covs), Y=recD$recruits, C=clim_covs, 
                     parents1=recD$parents1, parents2=recD$parents2,
                     G=G, gid=groups, tau=tau.set[j],
                     npreds=npreds, y_holdout=y_holdout, C_out=clim_holds,
                     parents1_out=parents1_out, parents2_out=parents2_out,
                     gid_out=groups_out, Gout=Gout)
    inits=list()
    inits[[1]]=list(a=rep(4,Yrs), a_mu=1, sig_a=1,
                    gint=rep(0,G), sig_G=1, u=0.4,
                    dd=-1,theta=1, b2=rep(0,length(clim_covs)))
    pars <- c("b2", "lppd")
    fitted <- stan(fit=mcmc_samples, data=datalist, pars=pars,
                   chains=1, iter=iterations, warmup=n.mcmc, 
                   init = list(inits[[1]]))
    long <- ggs(fitted)
    postmeans <- ddply(long, .(Parameter), summarise,
                       postmean = mean(value))
    tmpbetas <- postmeans[grep("b2", postmeans$Parameter), "postmean"]
    beta_cvs[s,] <- tmpbetas
    lppd[s] <- sum(postmeans[grep("lppd", postmeans$Parameter), "postmean"])
    if(is.na(lppd[s])){stop("bad lppd")}
  } # end C-V loop
  
  beta_regs[,,j] <- beta_cvs
  lppd2[,j] <- lppd

} # end regularization loop

## Sum across C-Vs
temp1 <- colSums(lppd2)
plot(log(1/tau.set), temp1, type="l")
temp2 <- colMeans(beta_regs)
matplot(t(temp2), type="l", lwd=2)
abline(h=0, lty=2)

# 
# 
 
# 
#   
#   
#   
# rng_seed <- 123
# sflist <-
#   mclapply(1:3, mc.cores=3,
#            function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
#                             seed=rng_seed, chains=1, chain_id=i, refresh=-1,
#                             iter=200, warmup=100, init=list(inits[[i]])))
# fit <- sflist2stanfit(sflist)
# saveRDS(fit, "recruitment_stanmcmc_PASM.RDS")


#     for(i in 1:max(long$Iteration)){
#       sublong <- subset(long, Iteration==i)
#       theta <- unlist(sublong[grep("theta", sublong$Parameter),"value"])
#       intercept <- unlist(sublong[grep("a_mu", sublong$Parameter),"value"])
#       group_eff <- unlist(sublong[grep("gint.", sublong$Parameter),"value"])
#       clim_eff <- unlist(sublong[grep("b2", sublong$Parameter),"value"])
#       crowd <- unlist(sublong[grep("dd", sublong$Parameter),"value"])
#       u <- sublong[grep("u", sublong$Parameter), ]
#       u <- as.numeric(subset(u, Parameter=="u")[,"value"])
#       trueP1 <- holdout$parents1*u + holdout$parents2*(1-u)
#       trueP2 <- sqrt(trueP1)
#       clim_outs <- holdout[,c(c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2"))]
#       climEff <- sum(as.numeric(clim_outs[1,])*clim_eff)
#       pred_df <- data.frame(obs=holdout$recruits, 
#                             p1=trueP1,
#                             p2=trueP2,
#                             intercept=rep(intercept, nrow(holdout)),
#                             crowd=rep(crowd, nrow(holdout)),
#                             climeffect=rep(climEff, nrow(holdout)))
#       mu <- with(pred_df, exp(intercept + crowd*p2 + climeffect))
#       lambda <- pred_df$p1*mu*theta
#       pred <- rnbinom(length(lambda), size = theta, mu = lambda)
#       # lppdsave[i] <- mean((pred-pred_df$obs)^2)
#       lppdsave[i] <- sum(dnbinom(x = pred_df$obs, size=theta, mu = lambda))
#     } # end MCMC prediction loop
