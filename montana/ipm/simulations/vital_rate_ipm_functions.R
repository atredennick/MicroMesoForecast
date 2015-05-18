####
#### Growth function
####
G <- function(v,u,W,Gpars,doYear,doSpp,climate){
  mu <- Gpars$intcpt[doSpp]+                    # intercept
        Gpars$intG[doSpp]+                      # group offset
        Gpars$slope[doSpp]*u+                   # size effect (slope)
        W*Gpars$nb[doSpp]+                      # crowding effect
        (W*Gpars$nbXsize[doSpp]*u)+                    # size-crowding interaction
        sum(Gpars$clim[,doSpp]*climate)+        # main climate effects, 
                                                #    including clim interactions
        sum(Gpars$slopeXclim[,doSpp]*climate[1:5]*u) # climate-size interaction effects

  sigma2 <- Gpars$sigma2.a[doSpp]*exp(Gpars$sigma2.b[doSpp]*mu)
  out <- dnorm(v,mu,sqrt(sigma2))
  out
}

####
#### Survival function: probability an individual of size 
#### u survives  (u is on log scale)
####
S <- function(u,W,Spars,doYear,doSpp,climate){
  mu <- Gpars$intcpt[doSpp]+                    # intercept
        Gpars$intG[doSpp]+                      # group offset
        Gpars$slope[doSpp]*u+                   # size effect (slope)
        W*Gpars$nb[doSpp]+                      # crowding effect
        (W*Gpars$nbXsize[doSpp]*u)+                    # size-crowding interaction
        sum(Gpars$clim[,doSpp]*climate)+        # main climate effects, 
                                                #    including clim interactions
        sum(Gpars$slopeXclim[,doSpp]*climate[1:5]*u) # climate-size interaction effects
  return(inv.logit(mu))
}


####
#### Number of recruits per area produced 
#### cover is stored in absolute area (cm^2)
####
get_rpa=function(Rpars,cover,climate){
  # cover is in m^2 per m^2; convert to % scale:
  cover2=cover*100
  # calculate recruits
  n_spp=length(cover)
  mu=rep(NA,n_spp)
  for(i in 1:n_spp){
    mu[i]=cover2[i]*exp(Rpars$intcpt[i]+sqrt(cover2[i])*Rpars$dd[i]+
                          Rpars$grpInt[i]+sum(Rpars$clim[,i]*climate)) 
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


