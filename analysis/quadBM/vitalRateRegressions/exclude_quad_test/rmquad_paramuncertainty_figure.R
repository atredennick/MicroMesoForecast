##  Script to make figure showing the effect of removing quads
##  on parameter uncertainty.

library(ggplot2)
library(ggmcmc)
library(reshape2)
library(plyr)
library(mgcv)

####
####  Read in model fits and summarize results
####
source("popgrowth_read_data.R")
spp <- unique(growD_all$Species)
num_quads <- numeric(length(spp))
for(i in 1:length(spp)){
  tmp <- subset(growD_all, Species==spp[i])
  num_quads[i] <- length(unique(tmp$quad))
}
num_quads <- as.data.frame(num_quads)
num_quads$species <- spp

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
####  Read in full fits
####
dir <- "../truncNormModel/"
files <- list.files(dir)
files <- files[grep(".RDS", files)]
all_alls <- data.frame(Parameter=NA, stddev=NA, id=NA, 
                       rmsquad=NA, rep=NA, species=NA)
for(i in 1:length(files)){
  tmp <- files[i]
  spp1 <- unlist(strsplit(tmp, "[_]"))[3]
  spp <- gsub(".RDS","",spp1)
  
  long <- readRDS(paste(dir,tmp, sep=""))
  short <- ddply(long, .(Parameter), summarise,
                 stddev = sd(value))
  climeff <- short[grep("b2", short$Parameter),]
  climeff$id <- substr(climeff$Parameter, 4, length(climeff$Parameter))
  climeff$id <- as.numeric(unlist(strsplit(climeff$id, split=']')))  
  climeff <- climeff[with(climeff, order(id)),]
  climeff$rmsquad <- 0
  climeff$rep <- 1
  climeff$species <- spp
  all_alls <- rbind(all_alls, climeff)
}
all_alls <- all_alls[2:nrow(all_alls),]

####
####  Aggregate results by species and rep
####
all_fits <- rbind(all_fits, all_alls)
agg_fits <- ddply(all_fits, .(Parameter, rmsquad, species, id), summarise,
                  avg_stddev = mean(stddev))
agg_fits2 <- agg_fits[which(agg_fits$id %in% c(1:7)),]
agg_fits2 <- merge(agg_fits2, num_quads, by="species")
agg_fits2$quadsfit <- with(agg_fits2, num_quads-as.numeric(rmsquad))


####
####  Fit model (need to do by species)
####
modBOGR <- gam(avg_stddev~quadsfit, data=subset(agg_fits2, species=="BOGR"),
               family=quasi(link="1/mu^2"))
modHECO <- gam(avg_stddev~quadsfit, data=subset(agg_fits2, species=="HECO"),
               family=quasi(link="1/mu^2"))
modPASM <- gam(avg_stddev~quadsfit, data=subset(agg_fits2, species=="PASM"),
               family=quasi(link="1/mu^2"))
modPOSE <- gam(avg_stddev~quadsfit, data=subset(agg_fits2, species=="POSE"),
               family=quasi(link="1/mu^2"))
summary(modBOGR)
summary(modHECO)
summary(modPASM)
summary(modPOSE)


####
#### Plot
####
ggplot(agg_fits2, aes(x=quadsfit, y=avg_stddev, color=Parameter))+
  geom_line(aes(group=Parameter), alpha=1, linetype=1, size=0.2)+
  geom_point(aes(group=Parameter), alpha=1, size=3)+
#   stat_smooth(method="loess", size=2, color="black", se=FALSE)+
#   geom_smooth(method="glm", formula=y~x, family=quasi(link="1/mu^2"),
#               size=2, color="black", se=FALSE)+
  facet_wrap("species", scales = "free", ncol=1)+
  theme_bw()+
  xlab("Number of quadrats fit")+
  ylab("Standard deviation of coefficient")

