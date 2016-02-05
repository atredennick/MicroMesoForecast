##
##  R Script to load and collate CV Regularization Results
##  from HPC user directory.
##
##  Author: Andrew Tredennick
##  Email: atredenn@gmail.com
##  Date created: 12-7-2015
##

rm(list=ls())

####
####  Load Libraries
####
library(plyr)
library(reshape2)
library(ggplot2)


####
####  Loop Through *.RDS Files and Load LPPD Results
####
do_spp <- "BOGR"
hpc_dir <- "/Volumes/A02046115/growth_oos/"
all_files <- list.files(hpc_dir)
lppd_files <- all_files[grep(".RDS", all_files)]
lppd_files <- lppd_files[grep(do_spp, lppd_files)]
lppd_mat <- matrix(NA, ncol=2, nrow=length(lppd_files))
for(i in 1:length(lppd_files)){
  tmp <- lppd_files[i]
  tmpsplit <- unlist(strsplit(tmp, "[.]"))
  tmp_grid_num <- unlist(strsplit(tmpsplit[1], "_"))[5]
  lppd_mat[i,1] <- as.numeric(tmp_grid_num)
  lppd_mat[i,2] <- as.numeric(readRDS(paste0(hpc_dir, tmp)))
}
lppd_mat <- as.data.frame(lppd_mat[order(lppd_mat[,1]),])
lppd_mat <- lppd_mat[complete.cases(lppd_mat),] 


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


##  Find out which grid ids didn't run on HPC; i.e., are missing
"%w/o%" <- function(x, y) x[!x %in% y] # x without y
allgrids <- c(1:n.grid)
rungrids <- lppd_mat[,1]
missedgrids <- allgrids %w/o% rungrids
paste0(missedgrids,collapse=",")
if(length(rungrids)+length(missedgrids)!=n.grid){ 
  stop("number of grids not matching up") }
if(length(missedgrids)>0){ 
  stop("missing some CV fits") }


####
####  Merge Grid and LPPD Results
####
colnames(lppd_mat) <- c("id", "lppd")
colnames(grid_df) <- c("sdev", "cv", "id")
lppd_mat2 <- merge(lppd_mat, grid_df, by="id")
lppd_cast <- dcast(lppd_mat2[,c("lppd", "sdev", "cv")], sdev~cv, value.var = "lppd")
score_cv_vec <- apply(lppd_cast,1,sum)
plot_df <- data.frame(betavar = sd_vec^2,
                      lppd = score_cv_vec)

opt_var <- plot_df[which(plot_df$lppd==max(plot_df$lppd)), "betavar"]

ggplot(data=plot_df, aes(x=betavar, y=lppd))+
  geom_line(size=1)+
  geom_vline(xintercept=opt_var, linetype=2)+
  xlab(bquote(sigma[beta]^2))+
  ylab("Log Predictive Score (lppd)")+
  theme_bw()


