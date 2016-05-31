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

data_path <- "./"
data_path <- "../../../../processed_data/"
allD <- readRDS(paste0(data_path,"rec_with_weather.RDS"))

year_ids <- sort(unique(allD$year))
yearD <- subset(allD, year!=year_ids[leave_out_year_id]) # remove leave-out-year


## Compile model outside of loop
recD <- subset(yearD, species==paste("R.",do_species,sep=""))


groups <- as.numeric(recD$group)
G <- length(unique(recD$group))
Yrs <- length(unique(recD$year))
yid <- as.numeric(as.factor(recD$year))

datalist <- list(N=nrow(recD), Yrs=Yrs, yid=yid,
                 Y=recD$recruits, 
                 parents1=recD$parents1, parents2=recD$parents2,
                 G=G, gid=groups)
pars=c("a_mu", "a", "u", "theta", "dd", "gint", "sig_a", "sig_G")
mcmc_samples <- stan(file="recruitment_noclimate.stan", data=datalist,pars=pars, chains=0)

inits=list()
inits[[1]]=list(a=rep(4,Yrs), a_mu=1, sig_a=1,
                gint=rep(0,G), sig_G=1, u=0.4,
                dd=-1,theta=1)
inits[[2]]=list(a=rep(0.5,Yrs), a_mu=0.2, sig_a=10,
                gint=rep(0,G), sig_G=0.1,  u=0.7,
                dd=-0.05,theta=1.5)
inits[[3]]=list(a=rep(1,Yrs), a_mu=0.5, sig_a=5,
                gint=rep(-0.1,G), sig_G=0.5,  u=0.5,
                dd=-1,theta=1)

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
saveRDS(long, paste("recruitment_stanmcmc_noclimate_", do_species, "_leaveout", year_ids[leave_out_year_id],".RDS", sep=""))

