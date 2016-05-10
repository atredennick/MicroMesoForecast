
library(reshape2)
library(plyr)
fitthin <- data.frame(Iteration=NA, Chain=NA, Parameter=NA,
                      value=NA, keep=NA, species=NA)
for(ispp in spp_list){
  fitlong <- readRDS(paste("../vitalRateRegs/validation/recruitment/fits_noclimate/recruitment_stanmcmc_noclimate_",ispp,"_leaveout",doYear, ".RDS", sep=""))
  fitlong$keep <- "no"
  keepseq <- seq(from = 1, to = nrow(fitlong), by = 10)
  fitlong[keepseq,"keep"] <- "yes"
  tmp <- subset(fitlong, keep=="yes")
  tmp$species <- ispp
  fitthin <- rbind(fitthin,tmp)
}
fitthin <- fitthin[2:nrow(fitthin),]

##  Break up MCMC into regression components
# Yearly intercepts
intercept_rec <- fitthin[grep("a", fitthin$Parameter),]
intercept_rec <- subset(intercept_rec, Parameter!="a_mu")
intercept_rec <- subset(intercept_rec, Parameter!="tau")
intercept_rec$yearid <- substr(intercept_rec$Parameter, 3, length(intercept_rec$Parameter))
intercept_rec$yearid <- unlist(strsplit(intercept_rec$yearid, split=']'))

# Mean intercept
interceptmu_rec <- fitthin[grep("a_mu", fitthin$Parameter),]

# Group effects
group_rec <- fitthin[grep("gint", fitthin$Parameter),]
group_rec$groupid <- substr(group_rec$Parameter, 6, length(group_rec$Parameter))
group_rec$groupid <- unlist(strsplit(group_rec$groupid, split=']'))

# Get parent density-dependent effects
densdep_rec <- fitthin[grep("dd", fitthin$Parameter),]

# Get mixing fraction u
mixfrac_rec <- fitthin[grep("u", fitthin$Parameter),]
mixfrac_rec <- subset(mixfrac_rec, Parameter!="a_mu")

# Get theta
theta_rec <- fitthin[grep("theta", fitthin$Parameter),]

# Get rid of big objects
rm(list = c("tmp","fitthin","fitlong"))

##  Define function to format survival coefficients
getRecCoefs <- function(doYear, groupnum){
  # Get random chain and iteration for this timestep
  tmp4chain <- subset(densdep_rec, species=="BOGR")
  randchain <- sample(x = tmp4chain$Chain, size = 1)
  randiter <- sample(x = tmp4chain$Iteration, size = 1)
  
  # Get random effects if doYear!=NA
  if(is.na(doYear)==FALSE){
    tmp_intercept <- subset(intercept_rec, yearid==doYear &
                              Iteration==randiter &
                              Chain==randchain)
  }
  
  # Set mean intercept and slope if doYear==NA
  if(is.na(doYear)==TRUE){
    tmp_intercept <- subset(interceptmu_rec, Iteration==randiter &
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
    tmp_group <- subset(group_rec, Iteration==randiter &
                          Chain==randchain &
                          groupid==groupnum)
    group_vec <- tmp_group$value
  }
  
  # Density dependence effect
  tmp_dd <- subset(densdep_rec, Iteration==randiter &
                            Chain==randchain)
  dd_vec <- tmp_dd$value
  
  # Mixing fraction
  tmp_u <- subset(mixfrac_rec, Iteration==randiter &
                           Chain==randchain)
  u_vec <- tmp_u$value
  
  # Theta
  tmp_theta <- subset(theta_rec, Iteration==randiter &
                             Chain==randchain)
  theta_vec <- tmp_theta$value
  
  ##  Collate all parameters for output
  Rpars=list(intcpt=intercept_vec, 
             grpInt=group_vec,
             dd=dd_vec,
             u=u_vec,
             theta=theta_vec)
  return(Rpars)
}

