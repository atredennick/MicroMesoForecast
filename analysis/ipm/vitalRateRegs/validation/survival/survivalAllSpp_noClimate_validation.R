##  The script takes 'grid_id' as a command line prompt to set spp and year left out
##  So, e.g., run as: "R CMD BATCH -1 survivalAllSpp_STAN.R" for species 1 and year 1 left out.

#clear everything, just to be safe 
rm(list=ls(all.names=TRUE))

##  Make species*year grid for calling model
sppList <- sort(c("BOGR","HECO","PASM","POSE"))
yearid <- 1:13
model_grid <- expand.grid(sppList, yearid)

## Set do_year and do_species for validation from command line prompt
args <- commandArgs(trailingOnly = F)
myargument <- args[length(args)]
myargument <- sub("-","",myargument)
grid_id <- as.numeric(myargument)
leave_out_year_id <- model_grid[grid_id, 2]
do_species <- as.character(model_grid[grid_id, 1])


#load libraries
library(rstan)
library(parallel)
library(reshape2)
library(ggmcmc)

# Read in data
data_path <- "./" # on HPC
data_path <- "../../../../processed_data/"
survD_all <- readRDS(paste0(data_path,"surv_with_weather.RDS"))

# Subset the data according to grid
year_ids <- sort(unique(survD_all$year))
yearD_all <- subset(survD_all, year!=year_ids[leave_out_year_id])

## Compile model 
survD <- subset(yearD_all, species==do_species)
groups <- as.numeric(survD$Group)
G <- length(unique(survD$Group))
Yrs <- length(unique(survD$year))
W <- cbind(survD$W, survD$W*log(survD$area))
yid <- as.numeric(as.factor(survD$year))

datalist <- list(N=nrow(survD), Yrs=Yrs, yid=yid,
                 Y=survD$survives, X=log(survD$area),
                 W=W, G=G, gid=groups)
pars=c("a_mu", "a", "b1_mu",  "b1",
       "w", "gint", "sig_a", "sig_b1", "sig_G")
mcmc_samples <- stan(file="survival_noclimate.stan", data=datalist, pars=pars, chains=0)

## Set reasonable initial values for three chains
inits <- list()
inits[[1]] <- list(a_mu=0, a=rep(0,Yrs), b1_mu=0.01, b1=rep(0.01,Yrs),
                   gint=rep(0,G), w=c(-0.05,-0.05), sig_b1=0.5, sig_a=0.5,
                   sig_G=0.5)
inits[[2]] <- list(a_mu=-0.1, a=rep(-0.1,Yrs), b1_mu=0.1, b1=rep(0.1,Yrs),
                   gint=rep(-0.1,G), w=c(-0.1,-0.1), sig_b1=0.2, sig_a=0.2,
                   sig_G=0.2)
inits[[3]] <- list(a_mu=0.05, a=rep(0.05,Yrs), b1_mu=0.05, b1=rep(0.05,Yrs),
                   gint=rep(0.1,G), w=c(-0.05,-0.05), sig_b1=0.3, sig_a=0.3,
                   sig_G=0.3)

##  Run MCMC in parallel
rng_seed <- 123
sflist <-
  mclapply(1:3, mc.cores=3,
          function(i) stan(fit=mcmc_samples, data=datalist, pars=pars,
                           seed=rng_seed, chains=1, chain_id=i, refresh=-1,
                           iter=2000, warmup=1000, init=list(inits[[i]])))
fit <- sflist2stanfit(sflist)
r_hats <- summary(fit)$summary[,10] 
write.csv(r_hats, paste("rhat_leaveout", year_ids[leave_out_year_id], "_noclimate_", do_species, ".csv", sep=""))
long <- ggs(fit)
saveRDS(long, paste("survival_stanmcmc_noclimate_", do_species, "_leaveout", year_ids[leave_out_year_id],".RDS", sep=""))
saveRDS(fit, outfile)
