model{
  #likelihood and process model
  for(i in 1:nObs){
    logit(gMu[i]) <- intYr[spp[i], yrs[i]] + intG[grp[i]] + beta[spp[i], yrs[i]]*X[i] + temp1[spp[i]]*TmeanSpr1[i] + temp2[spp[i]]*TmeanSpr2[i] + rain1[spp[i]]*ppt1[i] + rain2[spp[i]]*ppt2[i]  
    p[i] <- gMu[i] * tau[spp[i]]
    q[i] <- (1 - gMu[i]) * tau[spp[i]]
    C[i] ~ dbeta(p[i], q[i])
  }
  
  #priors
  for(s in 1:nSpp){
    betaSpp[s] ~ dnorm(0, 1e-6)
    temp1[s] ~ dnorm(temp1Mu, temp1Var)
    temp2[s] ~ dnorm(temp2Mu, temp2Var)
    rain1[s] ~ dnorm(rain1Mu, rain1Var)
    rain2[s] ~ dnorm(rain2Mu, rain2Var)
    intercept[s] ~ dnorm(0, 1e-6)
    intVaryY[s] ~ dgamma(0.001, 0.001)
    betaVar[s] ~ dgamma(0.001, 0.001)
    t0[s] ~ dnorm(0, .01) 
    tau[s] <- exp(t0[s])
    for(y in 1:nYrs){
      intYr[s,y] ~ dnorm(intercept[s], intVaryY[s])
      beta[s,y] ~ dnorm(betaSpp[s], betaVar[s])
    }
  }
  
  for(g in 1:nGrp){
    intG[g] ~ dnorm(0, intVarG)
  }
  
  temp1Mu ~ dnorm(0,1e-6)
  temp2Mu ~ dnorm(0,1e-6)
  rain1Mu ~ dnorm(0,1e-6)
  rain2Mu ~ dnorm(0,1e-6)
  temp1Var ~ dgamma(0.001, 0.001)
  temp2Var ~ dgamma(0.001, 0.001)
  rain1Var ~ dgamma(0.001, 0.001)
  rain2Var ~ dgamma(0.001, 0.001)
  intVarG ~ dgamma(2, 0.5) 
}