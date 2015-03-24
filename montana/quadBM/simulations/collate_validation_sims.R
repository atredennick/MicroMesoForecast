##  This script collates the leave-one-year-out validation simulations.

##  Author: Andrew Tredennick
##  Email:  atredenn@gmail.com

####
#### Load libraries ----------------------------------------
####
library('reshape2'); library('plyr'); library('ggplot2')


####
####  Read in simulation results 1-by-1 and combine --------
####
all_files <- list.files(path = "./validation_results/")
num_files <- length(all_files)
data_list <- list()
for(file in 1:num_files){
  now_file <- paste("./validation_results/", all_files[file], sep="")
  data_list[[file]] <- readRDS(now_file)
}

validation_sims <- melt(data_list, id.vars = c("quad", "species", "predicted_year", "sim")) 
#TODO: bring in observed data to calculate residuals

####
####  Quick plot -------------------------------------------
####
ggplot(validation_sims)+
  geom_boxplot(aes(x=as.character(predicted_year), y=value*100))+
  facet_wrap("species", scales = "free")