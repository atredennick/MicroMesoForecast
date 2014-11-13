model{
  #likelihood and process model
  for(i in 1:nObs){
    logit(p[i]) <- interceptMu + intG[grp[i]] + size[i]*betaMu 
    S[i] ~ dbern(p[i])
  }
  
  #priors
  for(g in 1:nGrp){
    intG[g] ~ dnorm(0, intVarG)
  }

  interceptMu ~ dnorm(0, 1e-6)
  betaMu ~ dnorm(0, 1e-6)
  intVarG ~ dgamma(2, 0.5) 
}