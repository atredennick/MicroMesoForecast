####
#### Growth function
####
G <- function(v,u,W,Gpars,doYear,doSpp,climate,Gscalers){
  clim_mains <- (climate - Gscalers[1:length(climate),"means"])/Gscalers[1:length(climate),"sds"]
  clim_size <- as.data.frame(as.matrix(apply(expand.grid(as.numeric(climate[1:5]),u), 1, prod)))
  clim_size$means <- rep(Gscalers[8:nrow(Gscalers),"means"], times=length(u))
  clim_size$sds <- rep(Gscalers[8:nrow(Gscalers),"sds"], times=length(u))
  clim_size$coefs <- rep(Gpars$slopeXclim[,doSpp], times=length(u))
  clim_size$climXsize <- with(clim_size, ((V1-means)/sds)*coefs)
  clim_size$midpoint <- rep(1:length(u), each=length(8:nrow(Gscalers)))
  climXsize <- ddply(clim_size, .(midpoint), summarise,
                     summed_effect = sum(climXsize))
  mu <- Gpars$intcpt[doSpp]+                    # intercept
        Gpars$intG[doSpp]+                      # group offset
        Gpars$slope[doSpp]*u+                   # size effect (slope)
        W*Gpars$nb[doSpp]+                      # crowding effect
        (W*Gpars$nbXsize[doSpp]*u)+                 # size-crowding interaction
        sum(Gpars$clim[,doSpp]*clim_mains)+        # main climate effects, 
                                                #    including clim interactions
        climXsize$summed_effect # climate-size interaction effects

  sigma2 <- Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  out <- dnorm(v,mu,sqrt(sigma2))
  out
}

####
#### Survival function: probability an individual of size 
#### u survives  (u is on log scale)
####
S <- function(u,W,Spars,doYear,doSpp,climate,Sscalers){
  clim_mains <- (climate - Sscalers[1:length(climate),"means"])/Sscalers[1:length(climate),"sds"]
  clim_size <- as.data.frame(as.matrix(apply(expand.grid(as.numeric(climate[1:5]),u), 1, prod)))
  clim_size$means <- rep(Sscalers[8:nrow(Sscalers),"means"], times=length(u))
  clim_size$sds <- rep(Sscalers[8:nrow(Sscalers),"sds"], times=length(u))
  clim_size$coefs <- rep(Gpars$slopeXclim[,doSpp], times=length(u))
  clim_size$climXsize <- with(clim_size, ((V1-means)/sds)*coefs)
  clim_size$midpoint <- rep(1:length(u), each=length(8:nrow(Sscalers)))
  climXsize <- ddply(clim_size, .(midpoint), summarise,
                     summed_effect = sum(climXsize))
  mu <- Spars$intcpt[doSpp]+                    # intercept
        Spars$intG[doSpp]+                      # group offset
        Spars$slope[doSpp]*u+                   # size effect (slope)
        W*Spars$nb[doSpp]+                      # crowding effect
        (W*Spars$nbXsize[doSpp]*u)  +                 # size-crowding interaction
        sum(Spars$clim[,doSpp]*clim_mains)+        # main climate effects, 
                                                #    including clim interactions
        climXsize$summed_effect # climate-size interaction effects
  return(inv.logit(mu))
}


####
#### Number of recruits per area produced 
#### cover is stored in absolute area (cm^2)
####
get_rpa=function(Rpars,cover,climate,Rscalers){
  # cover is in m^2 per m^2; convert to % scale:
  cover2=cover*100
  # calculate recruits
  n_spp=length(cover)
  mu=rep(NA,n_spp)
  for(i in 1:n_spp){
    clim_mains <- (climate - Rscalers[1:length(climate),"means"])/Rscalers[1:length(climate),"sds"]
    mu[i]=cover2[i]*exp(Rpars$intcpt[i]+sqrt(cover2[i])*Rpars$dd[i]+
                          Rpars$grpInt[i]+sum(Rpars$clim[i,]*clim_mains)) 
  }
  if(sum(is.na(mu))>0) browser() # stop for errors
  rpa=mu/(cover*Atotal)  # convert from number recruits to recruits per cm^2
  return(rpa)
}



####
#### Fecundity function: expected number of recruits of size y 
#### produced by a size x individual
#### The size distribution of recruits is on the log scale
####
f=function(v,u,Rpars,rpa,doSpp) { 
  nRecruits = rpa[doSpp]*exp(u)
  #probability of producing a seedling of size v
  tmp=dnorm(v,Rpars$sizeMean[doSpp],sqrt(Rpars$sizeVar[doSpp]))/(1-pnorm(-1.61,Rpars$sizeMean[doSpp],sqrt(Rpars$sizeVar[doSpp])))
  #number recruits of each size 
  f=nRecruits*tmp
  return(f)
}   


