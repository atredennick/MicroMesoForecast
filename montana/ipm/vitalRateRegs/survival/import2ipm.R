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

fitthin <- data.frame(Iteration=NA, Chain=NA, Parameter=NA,
                      value=NA, keep=NA, species=NA)
for(ispp in spp_list){
  fitlong <- readRDS(paste("../vitalRateRegs/survival/survival_stanmcmc_", ispp, ".RDS", sep=""))
  fitlong$keep <- "no"
  keepseq <- seq(from = 1, to = nrow(fitlong), by = 10)
  fitlong[keepseq,"keep"] <- "yes"
  tmp <- subset(fitlong, keep=="yes")
  tmp$species <- ispp
  fitthin <- rbind(fitthin,tmp)
}
fitthin <- fitthin[2:nrow(fitthin),]

##  Break up MCMC into regression components
# Climate effects
climeff <- fitthin[grep("b2", fitthin$Parameter),]

# Yearly cover (size) effects
coveff <- fitthin[grep(glob2rx("b1[*]"), fitthin$Parameter),]
coveff$yearid <- substr(coveff$Parameter, 4, length(coveff$Parameter))
coveff$yearid <- unlist(strsplit(coveff$yearid, split=']'))

# Mean cover effect
covermu <- fitthin[grep("b1_mu", fitthin$Parameter),]

# Yearly intercepts
intercept <- fitthin[grep("a", fitthin$Parameter),]
intercept <- subset(intercept, Parameter!="a_mu")
intercept <- subset(intercept, Parameter!="tau")
intercept$yearid <- substr(intercept$Parameter, 3, length(intercept$Parameter))
intercept$yearid <- unlist(strsplit(intercept$yearid, split=']'))

# Mean intercept
interceptmu <- fitthin[grep("a_mu", fitthin$Parameter),]

# Crowding effects
crowd <- fitthin[grep("w", fitthin$Parameter),]

# Group effects
group <- fitthin[grep("gint", fitthin$Parameter),]
group$groupid <- substr(group$Parameter, 6, length(group$Parameter))
group$groupid <- unlist(strsplit(group$groupid, split=']'))

## Get rid of big objects
rm(list = c("tmp","fitthin","fitlong"))

##  Define function to format survival coefficients
getSurvCoefs <- function(doYear, groupnum){
  # Get random chain and iteration for this timestep
  tmp4chain <- subset(climeff, species=="BOGR")
  randchain <- sample(x = tmp4chain$Chain, size = 1)
  randiter <- sample(x = tmp4chain$Iteration, size = 1)
  
  # Get random effects if doYear!=NA
  if(is.na(doYear)==FALSE){
    tmp_intercept <- subset(intercept, yearid==doYear &
                                       Iteration==randiter &
                                       Chain==randchain)
    tmp_size <- subset(coveff, yearid==doYear &
                               Iteration==randiter &
                               Chain==randchain)
  }
  
  # Set mean intercept and slope if doYear==NA
  if(is.na(doYear)==TRUE){
    tmp_intercept <- subset(interceptmu, Iteration==randiter &
                                         Chain==randchain)
    tmp_size <- subset(covermu, Iteration==randiter &
                                Chain==randchain)
  }
  size_vec <- tmp_size$value
  intercept_vec <- tmp_intercept$value
  
  ##  Now do group, climate, and competition fixed effects
  # Group effects
  if(is.na(groupnum)==TRUE){
    tmp_group <- rep(0,length(spp_list))
  }
  if(is.na(groupnum)==FALSE){
    tmp_group <- subset(group, Iteration==randiter &
                               Chain==randchain &
                               groupid==groupnum)
  }
  group_vec <- tmp_group$value
  
  # Climate effects
  tmp_clim <- subset(climeff, Iteration==randiter &
                              Chain==randchain)
  clim_mat <- matrix(tmp_clim$value, length(unique(tmp_clim$Parameter)), length(spp_list))
  
  # Crowding effect
  tmp_crowd <- subset(crowd, Iteration==randiter &
                             Chain==randchain)
  crowd_mat <- matrix(tmp_crowd$value, length(unique(tmp_crowd$Parameter)), length(spp_list))
  
  ##  Collate all parameters for output
  Spars=list(intcpt=intercept_vec, 
             intG=group_vec,
             slope=size_vec,
             nb=crowd_mat,
             clim=clim_mat)
  return(Spars)
}
