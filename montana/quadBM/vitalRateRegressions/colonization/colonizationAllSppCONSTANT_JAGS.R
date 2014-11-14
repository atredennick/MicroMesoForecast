model{
  #likelihood and process model
  for(i in 1:nObs){
    logit(p[i]) <- intercept[spp[i]] + intG[spp[i],grp[i]]
    C[i] ~ dbern(p[i])
  }
  
  #priors
  for(s in 1:nSpp){
    intercept[s] ~ dnorm(0, 1e-6)
    intVarG[s] ~ dgamma(2, 0.5) 
    for(g in 1:nGrp){
      intG[s,g] ~ dnorm(0, intVarG[s])
    }
  }
}