
library(reshape2)
library(plyr)
fitthin <- data.frame(Iteration=NA, Chain=NA, Parameter=NA,
                      value=NA, keep=NA, species=NA)
for(ispp in spp_list){
  fitlong <- readRDS(paste("../vitalRateRegs/recruitment/recruitment_stanmcmc_", ispp, ".RDS", sep=""))
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

# Yearly intercepts
intercept <- fitthin[grep("a", fitthin$Parameter),]
intercept <- subset(intercept, Parameter!="a_mu")
intercept <- subset(intercept, Parameter!="tau")
intercept$yearid <- substr(intercept$Parameter, 3, length(intercept$Parameter))
intercept$yearid <- unlist(strsplit(intercept$yearid, split=']'))

# Mean intercept
interceptmu <- fitthin[grep("a_mu", fitthin$Parameter),]

# Group effects
group <- fitthin[grep("gint", fitthin$Parameter),]
group$groupid <- substr(group$Parameter, 6, length(group$Parameter))
group$groupid <- unlist(strsplit(group$groupid, split=']'))

# Get parent density-dependent effects
densdep <- fitthin[grep("dd", fitthin$Parameter),]

# Get mixing fraction u
mixfrac <- fitthin[grep("u", fitthin$Parameter),]
mixfrac <- subset(mixfrac, Parameter!="a_mu")

# Get theta
theta <- fitthin[grep("theta", fitthin$Parameter),]

# Get rid of big objects
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
  }
  
  # Set mean intercept and slope if doYear==NA
  if(is.na(doYear)==TRUE){
    tmp_intercept <- subset(interceptmu, Iteration==randiter &
                              Chain==randchain)
  }
  intercept_vec <- tmp_intercept$value
  
  ##  Now do group, climate, and competition fixed effects
  # Group effects
  if(is.na(groupnum)==TRUE){
    tmp_group <- rep(0,length(spp_list))
    group_vec <- tmp_group
  }
  if(is.na(groupnum)==FALSE){
    tmp_group <- subset(group, Iteration==randiter &
                          Chain==randchain &
                          groupid==groupnum)
    group_vec <- tmp_group$value
  }
  
  # Climate effects
  tmp_clim <- subset(climeff, Iteration==randiter &
                       Chain==randchain)
  clim_mat <- matrix(tmp_clim$value, length(unique(tmp_clim$Parameter)), length(spp_list))
  
  # Density dependence effect
  tmp_dd <- subset(densdep, Iteration==randiter &
                            Chain==randchain)
  dd_vec <- tmp_dd$value
  
  # Mixing fraction
  tmp_u <- subset(mixfrac, Iteration==randiter &
                           Chain==randchain)
  u_vec <- tmp_u$value
  
  # Theta
  tmp_theta <- subset(theta, Iteration==randiter &
                             Chain==randchain)
  theta_vec <- tmp_theta$value
  
  ##  Collate all parameters for output
  Rpars=list(intcpt=intercept_vec, 
             grpInt=group_vec,
             dd=dd_vec,
             clim=clim_mat,
             u=u_vec,
             theta=theta_vec)
  return(Rpars)
}

