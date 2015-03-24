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
observed_data <- read.csv("../../speciesData/quadAllCover.csv")
combined_set <- merge(validation_sims, observed_data,
                      by.x=c("quad", "predicted_year", "species"), 
                      by.y=c("quad", "year", "Species"))
combined_set$residuals <- with(combined_set, value*100- propCover*100)

####
####  Quick plot -------------------------------------------
####
ggplot(combined_set)+
  geom_hline(aes(yintercept=0))+
  geom_boxplot(aes(x=as.character(predicted_year), y=residuals))+
  facet_wrap("species", scales = "free")+
  ylab("Residuals (predicted - observed)")


####
####  Write output as one file ----------------------------
####
out_set <- combined_set[,c("quad", "predicted_year", "species", "sim", "value", "propCover", "residuals")]
colnames(out_set) <- c("quad", "predicted_year", "species", "sim", "predicted_cover", "observed_cover", "residuals")
saveRDS(out_set, "qbm_loyo_validation_results.rds")



