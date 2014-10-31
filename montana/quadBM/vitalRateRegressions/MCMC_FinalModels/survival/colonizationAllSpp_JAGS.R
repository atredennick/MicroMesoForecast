model{
  #likelihood and process model
  for(i in 1:nObs){
    logit(p[i]) <- interceptMu 
    C[i] ~ dbern(p[i])
  }
  
#   #priors
#   for(g in 1:nGrp){
#     intG[g] ~ dnorm(0, intVarG)
#   }
#   
#   for(y in 1:nYrs){
#     intYr[y] ~ dnorm(interceptMu, intVarY)
#   }
  
  temp1 ~ dnorm(0,1e-6)
  temp2 ~ dnorm(0,1e-6)
  rain1 ~ dnorm(0,1e-6)
  rain2 ~ dnorm(0,1e-6)
  interceptMu ~ dnorm(0, 0.0001)
  intVarY ~ dgamma(0.001, 0.001)
  intVarG ~ dgamma(2, 0.5) 
}