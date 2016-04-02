##  The script takes 'do_species' as a command line prompt. So, e.g.,
##    run as: "R CMD BATCH -1 survivalAllSpp_STAN.R" for species 1.

#clear everything, just to be safe 
rm(list=ls(all=TRUE))

#load libraries
library(rstan)
library(parallel)

##  Read in lppd scores and selected prior climate stddevs
priors_df <- read.csv("../../../all_maxlppds.csv")
priors <- subset(priors_df, vital=="survival")

## Set do_year for validation from command line prompt
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
do_species <- as.numeric(myargument)
sppList=sort(c("BOGR","HECO","PASM","POSE"))
do_species <- 1 # for testing; remove for HPC runs

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

data_path <- "../../../speciesData/" #on local machine
# data_path <- "speciesData/" #on HPC server

for(spp in 1:length(sppList)){
  doSpp <- sppList[spp]
  
  if(doSpp == "BOGR"){
    sppD <- read.csv(paste(data_path, doSpp, "/edited/survD.csv", sep=""))
    sppD$species <- doSpp 
  }else{
    sppD <- read.csv(paste(data_path, doSpp, "/survD.csv", sep=""))
    sppD$species <- doSpp 
  }
  outD <- rbind(outD, sppD)
}

survD <- outD[2:nrow(outD),]

climD <- read.csv("../../../weather/Climate.csv") #on local machine
clim_vars <- c("pptLag", "ppt1", "ppt2", "TmeanSpr1", "TmeanSpr2")
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
colnames(crowd) <- c("W", "X", "species")

# Merge crowding and growth data
survD_all <- merge(survD, crowd, by=c("species", "X"))

## Compile model 
survD <- subset(survD_all, species==sppList[do_species]) #  grab a species

##  Create and scale interaction covariates
survD$ppt1TmeanSpr1 <- survD$ppt1*survD$TmeanSpr1
survD$ppt2TmeanSpr2 <- survD$ppt2*survD$TmeanSpr2
clim_vars_all <- c(clim_vars, "ppt1TmeanSpr1", "ppt2TmeanSpr2")
clim_covs <- survD[,clim_vars_all]
clim_covs <- scale(clim_covs, center = TRUE, scale = TRUE) # center and scale

##  Create objects for other effects
groups <- as.numeric(survD$Group)
G <- length(unique(survD$Group))
nyrs <- length(unique(survD$year))
W <- cbind(survD$W, survD$W*log(survD$area))
yid <- as.numeric(as.factor(survD$year))

# Get climate prior std dev
prior_stddev <- as.numeric(subset(priors, species==sppList[do_species])["prior_stdev"])

datalist <- list(N=nrow(survD), Yrs=nyrs, yid=yid,
                 Covs=ncol(clim_covs), Y=survD$survives, X=log(survD$area),
                 C=clim_covs, W=W, G=G, gid=groups, beta_tau=prior_stddev)
pars=c("a_mu", "a", "b1_mu",  "b1", "b2",
       "w", "gint")

mcmc_samples <- stan(file="survival.stan", data=datalist, pars=pars, chains=0)

## Set reasonable initial values for three chains
inits <- list()
inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                   gint=rep(0,G), w=c(-0.05,-0.05), sig_b1=0.5, sig_a=0.5,
                   sig_G=0.5, b2=rep(0,ncol(clim_covs)))
inits[[2]] <- list(a_mu=-0.1, a=rep(-0.1,Yrs), b1_mu=0.1, b1=rep(0.1,Yrs),
                   gint=rep(-0.1,G), w=c(-0.1,-0.1), sig_b1=0.2, sig_a=0.2,
                   sig_G=0.2, b2=rep(0.1,ncol(clim_covs)))
inits[[3]] <- list(a_mu=0.05, a=rep(0.05,Yrs), b1_mu=0.05, b1=rep(0.05,Yrs),
                   gint=rep(0.1,G), w=c(-0.05,-0.05), sig_b1=0.3, sig_a=0.3,
                   sig_G=0.3, b2=rep(0.05,ncol(clim_covs)))

##  Run MCMC in parallel
rng_seed <- 123
sflist <-
  mclapply(1:3, mc.cores=3,
          function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                           seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                           iter=2000, warmup=1000, init=list(inits[[i]])))
fit <- sflist2stanfit(sflist)
longfit <- ggs(fit) # convert StanFit --> data frame
outfile <- paste("survival_stanmcmc_", sppList[do_species], ".RDS", sep="")
saveRDS(fit, outfile)
