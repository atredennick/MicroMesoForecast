model{
  #likelihood and process models
  for(i in 1:nObs){
    Y[i] ~ dnegbin(q[i], theta[spp[i]])
    q[i] <- theta[spp[i]] / (theta[spp[i]]+lambda[i])
    lambda[i] <- trueP1[i]*fecundity[i]
    log(fecundity[i]) <- intcpt.mu[spp[i]] + intG[grp[i], spp[i]] + dd[spp[i]]*trueP2[i]
    trueP1[i] <- parents1[i]*u[spp[i]] + parents2[i]*(1-u[spp[i]])
    trueP2[i] <- sqrt(trueP1[i])
  }
  
  #priors
  for(s in 1:nSpp){
    u[s] ~ dunif(0,1)
    theta[s] ~ dgamma(0.001, 0.001)
    g.tau[s] ~ dgamma(2, 0.5)
    intcpt.mu[s] ~ dnorm(0, 1e-6)
    intcpt.tau[s] ~ dgamma(0.001, 0.001)
    dd[s] ~ dnorm(0, 1e-6)
    for(g in 1:nGrp){
      intG[g,s] ~ dnorm(0, g.tau[s])
    }#end group loop 
  }#end species loop 
}#end model
