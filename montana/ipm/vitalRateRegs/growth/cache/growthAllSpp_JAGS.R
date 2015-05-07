model{
  #process model and likelihood
  for(i in 1:nObs){
    mu[i] <- intYr[spp[i],yrs[i]] + intG[spp[i],grp[i]] + beta[spp[i],yrs[i]]*X[i] + 
             nb[spp[i]]*crowd[i] + temp1[spp[i]]*TmeanSpr1[i] + 
             temp2[spp[i]]*TmeanSpr2[i] + rain1[spp[i]]*ppt1[i] + 
             rain2[spp[i]]*ppt2[i] + rainlag[spp[i]]*pptlag[i]
    tau2[i] <- 1/(tau[spp[i]]*exp(tauSize[spp[i]]*mu[i])) 
    tau3[i] <- max(tau2[i],0.00000001)  
    Y[i] ~ dnorm(mu[i], tau3[spp[i]])
  }
  
  #priors
  for(s in 1:nSpp){
#     tau[s] ~ dgamma(0.5, 0.5)
    tau[s] ~ dnorm(0,0.001)
    tauSize[s] ~ dnorm(0,0.001)
    betaSpp[s] ~ dnorm(0,0.001)
    nb[s] ~ dunif(-5, 5)
    temp1[s] ~ dunif(-5, 5)
    temp2[s] ~ dunif(-5, 5)
    rain1[s] ~ dunif(-5, 5)
    rain2[s] ~ dunif(-5, 5)
    rainlag[s] ~ dunif(-5, 5)
    intercept[s] ~ dnorm(0, 0.001)
    intVaryY[s] ~ dgamma(0.5, 0.5)
    betaVar[s] ~ dgamma(0.5, 0.5)
    intVarG[s] ~ dgamma(0.001, 0.001) 
    for(y in 1:nYrs){
      intYr[s,y] ~ dnorm(intercept[s], intVaryY[s])
      beta[s,y] ~ dnorm(betaSpp[s], betaVar[s])
    }
    for(g in 1:nGrp){
      intG[s,g] ~ dnorm(0, intVarG[s])
    }
  }
  
#   temp1Mu ~ dnorm(0,1e-6)
#   temp2Mu ~ dnorm(0,1e-6)
#   rain1Mu ~ dnorm(0,1e-6)
#   rain2Mu ~ dnorm(0,1e-6)
#   rainlagMu ~ dnorm(0,1e-6)
#   temp1Var ~ dgamma(0.001, 0.001)
#   temp2Var ~ dgamma(0.001, 0.001)
#   rain1Var ~ dgamma(0.001, 0.001)
#   rain2Var ~ dgamma(0.001, 0.001)
#   rainlagVar ~ dgamma(0.001, 0.001)
}