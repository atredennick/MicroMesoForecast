##  Script to source all PASM regressions
setwd("growth/")
source("growthPASM_STAN.R")
setwd("../survival/")
source("survivalPASM_STAN.R")
setwd("../recruitment/")
source("recruitmentPASM_STAN.R")