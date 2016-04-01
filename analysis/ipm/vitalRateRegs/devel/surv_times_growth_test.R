##  Script to test sensitivity of predicted cover with and without
##  climate by size interaction terms.

rm(list=ls())

# Libraries
library(lme4)

# Set species
do_species <- "HECO"

# Source growth regressions (1=without climXsize, 2=with)
source("growthAllSpp_GLM.R")

# Source survival regressions (1=without climXsize, 2=with)
source("survivalAllSpp_GLM.R")

# Predict a quadrat-year
reps <- 500
comp1 <- numeric(reps)
comp1[] <- NA
comp2 <- comp1
for(i in 1:reps){
  quads <- unique(survD$quad)
  years <- unique(survD$year)
  do_year <- sample(years,1)
  do_quad <- sample(quads,1)
  pred_data <- subset(survD, quad==do_quad & year==do_year)
  growth_data <- subset(growD, quad==do_quad & year==do_year)
  if( nrow(pred_data) > 0 ){
    # S1 <- rbinom(nrow(pred_data), size = 1, prob = predict(surv_out1, pred_data, type="response"))
    S1 <- predict(surv_out1, pred_data, type="response")
    S2 <- predict(surv_out2, pred_data, type="response")
    pred_datag <- pred_data
    colnames(pred_datag)[which(colnames(pred_datag)=="area")] <- "area.t0"
    G1 <- exp(predict(growth_out1, pred_datag, type="response"))
    G2 <- exp(predict(growth_out2, pred_datag, type="response"))
    pred_cover1 <- sum(S1*G1)
    pred_cover2 <- sum(S2*G2)
    obs_cover <- sum(growth_data$area.t1)
    comp1[i] <- (pred_cover1-obs_cover)/100
    comp2[i] <- (pred_cover2-obs_cover)/100
  }
}
par(mfrow=c(1,2))
hist(comp1, main="W/out sizeXclimate Interactions", xlab="Pred-Observed Cover (%)")
hist(comp2, main="With sizeXclimate Interactions", xlab="Pred-Observed Cover (%)")
