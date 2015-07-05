##  Script to make figure showing the effect of removing quads
##  on parameter uncertainty.

library(ggplot2)
library(ggmcmc)
library(reshape2)
library(plyr)

####
####  Read in model fits and summarize results
####

all_files <- list.files("./results/")
all_fits <- data.frame(Parameter=NA, stddev=NA, id=NA, 
                       rmsquad=NA, rep=NA, species=NA)
for(i in 1:length(all_files)){
  tmp <- all_files[i]
  spp <- substring(tmp, 1, 4)
  qsrm1 <- unlist(strsplit(tmp, "[_]"))[2]
  qsrm <- gsub("numout", "", qsrm1)
  rep1 <- unlist(strsplit(tmp, "[_]"))[3]
  rep2 <- gsub(".RDS","",rep1)
  rep <- gsub("rep", "", rep2)
  
  long <- readRDS(paste("./results/",tmp, sep=""))
  short <- ddply(long, .(Parameter), summarise,
                 stddev = sd(value))
  climeff <- short[grep("b2", short$Parameter),]
  climeff$id <- substr(climeff$Parameter, 4, length(climeff$Parameter))
  climeff$id <- as.numeric(unlist(strsplit(climeff$id, split=']')))  
  climeff <- climeff[with(climeff, order(id)),]
  climeff$rmsquad <- qsrm
  climeff$rep <- rep
  climeff$species <- spp
  all_fits <- rbind(all_fits, climeff)
} # end file loop

all_fits <- all_fits[2:nrow(all_fits),]


####
####  Aggregate results by species and rep
####

agg_fits <- ddply(all_fits, .(Parameter, rmsquad, species), summarise,
                  avg_stddev = mean(stddev))
ggplot(agg_fits, aes(x=as.numeric(rmsquad), y=avg_stddev))+
  geom_line(aes(group=Parameter), alpha=0.3)+
#   geom_point(aes(shape=Parameter), alpha=0.3, size=2)+
  stat_smooth(size=2, color="black")+
  facet_wrap("species", scales = "free")+
  theme_bw()
