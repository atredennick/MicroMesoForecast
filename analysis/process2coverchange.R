##  Script to calculate the percent cover change in a quad-year attributable
##  to growth, survival, or recruitment. Then we compare these to the mean
##  absolute error from the one-step validation results.

library(plyr)
library(reshape2)
library(ggplot2)

do_species <- "POSE"


####
####  Bring in observation data
####

##  Quadrat data
quad_data <- read.csv("./speciesData/quadAllCover.csv")
quad_data <- subset(quad_data, Species == do_species)

##  Growth data
file <- paste("./speciesData/", do_species, "/growDnoNA.csv", sep="")
if(do_species=="BOGR"){
  file <- paste("./speciesData/", do_species, "/edited/growDnoNA.csv", sep="")
}
grow_data <- read.csv(file)

##  Survival data  
file <- paste("./speciesData/", do_species, "/survD.csv", sep="")
if(do_species=="BOGR"){
  file <- paste("./speciesData/", do_species, "/edited/survD.csv", sep="")
}
surv_data <- read.csv(file)

##  Recruitment data
file <- paste("./speciesData/", do_species, "/recArea_quads_not_removed.csv", sep="")
if(do_species=="BOGR"){
  file <- paste("./speciesData/", do_species, "/edited/recArea_quads_not_removed.csv", sep="")
}
rec_data <- read.csv(file)


####
####  Calculate area changes from each process
####
grow_changes <- ddply(grow_data, .(quad, year), summarise,
                      grow_diff = sum(area.t1)-sum(area.t0))
surv_changes <- ddply(subset(surv_data, survives==0), .(quad, year), summarise,
                      mort_loss = -(sum(area)))
rec_changes <- ddply(rec_data, .(quad, year), summarise,
                     rec_gain = sum(recArea))
process_changes <- merge(grow_changes, surv_changes, all = TRUE)
process_changes <- merge(process_changes, rec_changes, all = TRUE)
process_changes[which(is.na(process_changes$grow_diff)==TRUE), "grow_diff"] <- 0
process_changes[which(is.na(process_changes$mort_loss)==TRUE), "mort_loss"] <- 0
process_changes[which(is.na(process_changes$rec_gain)==TRUE), "rec_gain"] <- 0
process_changes$totalAchange <- rowSums(process_changes[,c(3,4,5)])

####
####  Get total percent cover changes from quad data as check
####

##  Create lag cover variable
tmp <- quad_data[,c("quad","year","totCover")]
tmp$year <- tmp$year+1
names(tmp)[3]="lag.cover"
quad_data=merge(quad_data,tmp,all.x=T)
quad_data$totalAchange_quad <- quad_data$totCover - quad_data$lag.cover
###############################################################
## It compares well, so moving on just using individual data ##
###############################################################


####
####  Calculate contribution of each vital rate to changes
####  based on absolute magnitude of area change
####
process_changes$abs_tot_change <- rowSums(abs(process_changes[,c(3,4,5)]))
process_changes$growth_cont <- abs(process_changes$grow_diff) / process_changes$abs_tot_change
process_changes$surv_cont <- abs(process_changes$mort_loss) / process_changes$abs_tot_change
process_changes$rec_cont <- abs(process_changes$rec_gain) / process_changes$abs_tot_change


####
####  Plot
####

##  Make new quad-year variable
process_changes$quad_year <- paste(process_changes$quad, process_changes$year, sep="")

##  Bring in mean absolute error data from QBM forecasts
file <- paste("./quadBM/simulations/results/one_step_validation/", do_species, 
              "_sim_cover_1step_ahead_year.RDS", sep="")
onesteps <- readRDS(file)
onestep_agg <- ddply(onesteps, .(quad, year), summarise,
                     mae = mean(abs(cover.t1-cover.t0)),
                     sd = sd(cover.t1))
onestep_agg$year <- onestep_agg$year-1
plot_data <- merge(process_changes, onestep_agg, all = TRUE)

ggplot(plot_data, aes(x=growth_cont, y=mae))+
  geom_point()
