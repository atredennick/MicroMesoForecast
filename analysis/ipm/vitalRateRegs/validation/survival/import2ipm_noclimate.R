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

fitthin <- list()
for(ispp in spp_list){
  fitlong <- readRDS(paste("../vitalRateRegs/validation/survival/fits_noclimate/survival_stanmcmc_noclimate_",ispp,"_leaveout",doYear, ".RDS", sep=""))
  #   fitlong$keep <- "no"
  #   keepseq <- seq(from = 1, to = nrow(fitlong), by = 10)
  #   fitlong[keepseq,"keep"] <- "yes"
  #   tmp <- subset(fitlong, keep=="yes")
  #   tmp$species <- ispp
  fitlong$species <- ispp
  fitthin <- rbind(fitthin,fitlong)
}

##  Break up MCMC into regression components
# Yearly cover (size) effects
coveff_surv <- fitthin[grep(glob2rx("b1[*]"), fitthin$Parameter),]
coveff_surv$yearid <- substr(coveff_surv$Parameter, 4, length(coveff_surv$Parameter))
coveff_surv$yearid <- as.numeric(unlist(strsplit(coveff_surv$yearid, split=']')))

# Mean cover effect
covermu_surv <- fitthin[grep("b1_mu", fitthin$Parameter),]
coversig_surv <- fitthin[grep("sig_b1", fitthin$Parameter),]

# Yearly intercepts
intercept_surv <- fitthin[grep("a", fitthin$Parameter),]
intercept_surv <- subset(intercept_surv, Parameter!="a_mu")
intercept_surv <- subset(intercept_surv, Parameter!="tau")
intercept_surv$yearid <- substr(intercept_surv$Parameter, 3, length(intercept_surv$Parameter))
intercept_surv$yearid <- unlist(strsplit(intercept_surv$yearid, split=']'))

# Mean intercept
interceptmu_surv <- fitthin[grep("a_mu", fitthin$Parameter),]
interceptsig_surv <- fitthin[grep("sig_a", fitthin$Parameter),]

# Crowding effects
crowd_surv <- fitthin[grep("w", fitthin$Parameter),]

# Group effects
group_surv <- fitthin[grep("gint", fitthin$Parameter),]
group_surv$groupid <- substr(group_surv$Parameter, 6, length(group_surv$Parameter))
group_surv$groupid <- unlist(strsplit(group_surv$groupid, split=']'))

## Get rid of big objects
rm(list = c("tmp","fitthin","fitlong"))

##  Define function to format survival coefficients
getSurvCoefs <- function(doYear, groupnum){
  # Get random chain and iteration for this timestep
  tmp4chain <- subset(coveff_surv, species=="BOGR")
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
    
    size_vec <- tmp_size$value
    intercept_vec <- tmp_intercept$value
  }
  
  # Set mean intercept and slope if doYear==NA
  if(is.na(doYear)==TRUE){
    tmp_intercept_mu <- subset(interceptmu_surv, Iteration==randiter &
                                 Chain==randchain)
    tmp_intercept_sig <- subset(interceptsig_surv, Iteration==randiter &
                                  Chain==randchain)
    tmp_intercept <- rnorm(4, tmp_intercept_mu$value, tmp_intercept_sig$value)
    tmp_size_mu <- subset(covermu_surv, Iteration==randiter &
                            Chain==randchain)
    tmp_size_sig <- subset(covermu_surv, Iteration==randiter &
                             Chain==randchain)
    tmp_size <- rnorm(4, tmp_size_mu$value, tmp_size_sig$value)
    
    size_vec <- tmp_size
    intercept_vec <- tmp_intercept
  }
  
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
  
  # Crowding effect
  tmp_crowd <- subset(crowd_surv, Iteration==randiter &
                             Chain==randchain)
  crowd_mat <- matrix(tmp_crowd$value, length(unique(tmp_crowd$Parameter)), length(spp_list))
  
  ##  Collate all parameters for output
  Spars=list(intcpt=intercept_vec, 
             intG=group_vec,
             slope=size_vec,
             nb=crowd_mat[1,],
             nbXsize=crowd_mat[2,])
  return(Spars)
}
