# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Test script to pull in growth parameters for IPM #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#Several arguments are required
#' @param doYear  Specific climate year to pull year random effects (19xx)
#' @param mcDraw  A numeric scalar or vector for the row(s) of MCMC to draw parameters from 
#' @param doSpp   A character scalar for the current speicies
#' @param group   A numeric scalar with group indicator

library(reshape2)
library(plyr)
library(ggmcmc)

fitthin <- data.frame(Iteration=NA, Chain=NA, Parameter=NA,
                      value=NA, keep=NA, species=NA)
for(ispp in spp_list){
  fitlong <- ggs(readRDS(paste("../vitalRateRegs/survival/survival_stanmcmc_", ispp, ".RDS", sep="")))
  fitlong$keep <- "no"
  keepseq <- seq(from = 1, to = nrow(fitlong), by = 1)
  fitlong[keepseq,"keep"] <- "yes"
  tmp <- subset(fitlong, keep=="yes")
  tmp$species <- ispp
  fitthin <- rbind(fitthin,tmp)
}
fitthin <- fitthin[2:nrow(fitthin),]

##  Break up MCMC into regression components
# Climate effects
climeff_surv <- fitthin[grep("b2", fitthin$Parameter),]

# Yearly cover (size) effects
coveff_surv <- fitthin[grep("b1[.]", fitthin$Parameter),]
coveff_surv$yearid <- substr(coveff_surv$Parameter, 4, length(coveff_surv$Parameter))
coveff_surv$yearid <- unlist(strsplit(coveff_surv$yearid, split='[.]'))

# Mean cover effect
covermu_surv <- fitthin[grep("b1_mu", fitthin$Parameter),]

# Yearly intercepts
intercept_surv <- fitthin[grep("a", fitthin$Parameter),]
intercept_surv <- subset(intercept_surv, Parameter!="a_mu")
intercept_surv <- subset(intercept_surv, Parameter!="tau")
intercept_surv$yearid <- substr(intercept_surv$Parameter, 3, length(intercept_surv$Parameter))
intercept_surv$yearid <- unlist(strsplit(intercept_surv$yearid, split='[.]'))

# Mean intercept
interceptmu_surv <- fitthin[grep("a_mu", fitthin$Parameter),]

# Crowding effects
crowd_surv <- fitthin[grep("w", fitthin$Parameter),]

# Group effects
group_surv <- fitthin[grep("gint", fitthin$Parameter),]
group_surv$groupid <- substr(group_surv$Parameter, 6, length(group_surv$Parameter))
group_surv$groupid <- unlist(strsplit(group_surv$groupid, split='[.]'))

## Get rid of big objects
rm(list = c("tmp","fitthin","fitlong"))

##  Define function to format survival coefficients
getSurvCoefs <- function(doYear, groupnum){
  # Get random chain and iteration for this timestep
  tmp4chain <- subset(climeff_surv, species=="BOGR")
  randchain <- sample(x = tmp4chain$Chain, size = 1)
  randiter <- sample(x = tmp4chain$Iteration, size = 1)
  
  # Get random effects if doYear!=NA
  if(is.na(doYear)==FALSE){
    tmp_intercept <- subset(intercept_surv, yearid==doYear &
                                       Iteration==randiter &
                                       Chain==randchain)
    tmp_size <- subset(coveff_surv, yearid==doYear &
                               Iteration==randiter &
                               Chain==randchain)
  }
  
  # Set mean intercept and slope if doYear==NA
  if(is.na(doYear)==TRUE){
    tmp_intercept <- subset(interceptmu_surv, Iteration==randiter &
                                         Chain==randchain)
    tmp_size <- subset(covermu_surv, Iteration==randiter &
                                Chain==randchain)
  }
  size_vec <- tmp_size$value
  intercept_vec <- tmp_intercept$value
  
  ##  Now do group, climate, and competition fixed effects
  # Group effects
  if(is.na(groupnum)==TRUE){
    tmp_group <- rep(0,length(spp_list))
    group_vec <- tmp_group
  }
  if(is.na(groupnum)==FALSE){
    tmp_group <- subset(group_surv, Iteration==randiter &
                               Chain==randchain &
                               groupid==groupnum)
    group_vec <- tmp_group$value
  }
  
  # Climate effects
  tmp_clim <- subset(climeff_surv, Iteration==randiter &
                              Chain==randchain)
  clim_mat <- matrix(tmp_clim$value, length(unique(tmp_clim$Parameter)), length(spp_list))
  
  # Crowding effect
  tmp_crowd <- subset(crowd_surv, Iteration==randiter &
                             Chain==randchain)
  crowd_mat <- matrix(tmp_crowd$value, length(unique(tmp_crowd$Parameter)), length(spp_list))
  
  ##  Collate all parameters for output
  Spars=list(intcpt=intercept_vec, 
             intG=group_vec,
             slope=size_vec,
             nb=crowd_mat[1,],
             nbXsize=crowd_mat[2,],
             clim=clim_mat)
  return(Spars)
}
