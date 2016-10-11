##  Script to plot and save lppd results from ridge regression

rm(list=ls())


####
####  Load Libraries
####
library(plyr)
library(reshape2)
library(ggplot2)
library(ggthemes)



####
####  Recreate Sdev*CV Grid
####
n.beta <- 24
sd_vec <- seq(0.1,1.5,length.out = n.beta)
K <- 13 # number of years
cv.s2.grid <- expand.grid(1:n.beta,1:K)
n.grid <- dim(cv.s2.grid)[1]
fold.idx.mat <- matrix(1:13,ncol=K)
grid_df <- as.data.frame(cv.s2.grid)
grid_df$id <- rownames(grid_df)



####
####  Read in lppd Results
####
path2growth <- "/Volumes/A02046115/growth_oos/"
path2surv <- "/Volumes/A02046115/surv_oos/"
path2rec <- "/Volumes/A02046115/rec_oos/"
path2quad <- "/Volumes/A02046115/qbm_oos/"
allpaths <- c(path2growth, path2surv, path2rec, path2quad)
allvitals <- c("growth", "survival", "recruitment", "cover")
species <- c("BOGR", "HECO", "PASM", "POSE")
count <- 1

plot_df <- data.frame(betavar=NA, lppd=NA, species=NA, vital=NA)
for(do_path in allpaths){
  for(do_spp in species){
    all_files <- list.files(do_path)
    lppd_files <- all_files[grep(".RDS", all_files)]
    lppd_files <- lppd_files[grep(do_spp, lppd_files)]
    lppd_mat <- matrix(NA, ncol=2, nrow=length(lppd_files))
    for(i in 1:length(lppd_files)){
      tmp <- lppd_files[i]
      tmpsplit <- unlist(strsplit(tmp, "[.]"))
      tmp_grid_num <- unlist(strsplit(tmpsplit[1], "_"))[5]
      lppd_mat[i,1] <- as.numeric(tmp_grid_num)
      lppd_mat[i,2] <- as.numeric(readRDS(paste0(do_path, tmp)))
    }
    lppd_mat <- as.data.frame(lppd_mat[order(lppd_mat[,1]),])
    lppd_mat <- lppd_mat[complete.cases(lppd_mat),]  
    
    ####  Merge Grid and LPPD Results
    colnames(lppd_mat) <- c("id", "lppd")
    colnames(grid_df) <- c("sdev", "cv", "id")
    lppd_mat2 <- merge(lppd_mat, grid_df, by="id")
    lppd_cast <- dcast(lppd_mat2[,c("lppd", "sdev", "cv")], sdev~cv, value.var = "lppd")
    score_cv_vec <- apply(lppd_cast,1,sum)
    tmp_df <- data.frame(betavar = sd_vec^2,
                          lppd = score_cv_vec,
                          species = do_spp,
                          vital = allvitals[count])
    plot_df <- rbind(plot_df, tmp_df)
  } # end species loop
  count <- count+1 #for tracking vital rates
} # end path loop

plot_df <- plot_df[2:nrow(plot_df),]
max_df <- ddply(plot_df, .(species, vital), summarise,
                max_lppd = max(lppd))
max_df2 <- plot_df[which(plot_df$lppd%in%max_df$max_lppd),]

####
####  Make Plot
####
ggplot(plot_df, aes(x=betavar, y=lppd))+
  geom_line()+
  geom_point()+
  geom_vline(data=max_df2, aes(xintercept=betavar), color="red", linetype=2)+
  facet_wrap(vital~species, scales = "free")+
  xlab(bquote(sigma[beta]^2))+
  ylab("Log Predictive Score (lppd)")+
  theme_few()

ggsave("/Users/atredenn/Repos/MicroMesoForecast/manuscript/components/lppds_suppfig.png",width = 8.5, height=6, units="in", dpi=72)



####
####  Save Max LPPD Values and Variances
####
max_df2$prior_stdev <- round(sqrt(max_df2$betavar),2)
write.csv(max_df2, "/Users/atredenn/Repos/MicroMesoForecast/analysis/all_maxlppds.csv")
