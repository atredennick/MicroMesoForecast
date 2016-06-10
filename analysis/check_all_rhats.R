##  Script to plot all r-hats for vital rate model parameters

rm(list=ls(all.names = TRUE))

library(ggmcmc)

vitals <- c("survival", "growth", "recruitment")

####
####  VALIDATION FITS
####
bad_rhats <- list()
# Vital rate models
for(do_vital in vitals){
  indir <- paste0("./ipm/vitalRateRegs/validation/", do_vital,"/rhats/")
  all_files <- list.files(indir)
  for(i in 1:length(all_files)){
    do_file <- read.csv(paste0(indir,all_files[i]))
    do_file <- subset(do_file, X!="lp__")
    if(max(do_file$x) > 1.1){
      tmp_df <- data.frame(file=paste0(indir,all_files[i]), max_rhat=max(do_file$x), parameter=do_file[which(do_file$x==max(do_file$x)),"X"])
      bad_rhats <- rbind(bad_rhats, tmp_df)
      cat("BAD: rhat > 1.1 \n")
    } # end max 1.1 if/then
    if(max(do_file$x) < 1.1){
      cat("GOOD: rhat < 1.1 \n")
    } # end max 1.1 if/then
  }
}


# Percent cover models
indir <- "./quadBM/vitalRateRegressions/truncNormModel/validation/rhats/"
all_files <- list.files(indir)
for(i in 1:length(all_files)){
  do_file <- read.csv(paste0(indir,all_files[i]))
  do_file <- subset(do_file, X!="lp__")
  if(max(do_file$x) > 1.1){
    tmp_df <- data.frame(file=paste0(indir,all_files[i]), max_rhat=max(do_file$x), parameter=do_file[which(do_file$x==max(do_file$x)),"X"])
    bad_rhats <- rbind(bad_rhats, tmp_df)
    cat("BAD: rhat > 1.1 \n")
  } # end max 1.1 if/then
  if(max(do_file$x) < 1.1){
    cat("GOOD: rhat < 1.1 \n")
  } # end max 1.1 if/then
}

bad_rhats



####
####  EXAMINE TRACEPLOTS OF HIGH R_HAT PARAMETERS
####
pdf(file =  "traceplots.pdf")
for(i in 1:nrow(bad_rhats)) {
  if(length(grep("noclimate", as.character(bad_rhats[i,"file"])))==0) {
    splits <- strsplit(as.character(bad_rhats[i,"file"]), "/")
    tmp_dir <- paste0(splits[[1]][1:5], collapse = "/")
    file_split1 <- strsplit(splits[[1]][7],"_")
    file_split2 <- strsplit(file_split1[[1]][3], "[.]")
    if(splits[[1]][2]=="ipm") {
      tmp_file <- paste0("recruitment_stanmcmc_", file_split2[[1]][1], "_", file_split1[[1]][2], ".RDS")
    }
    if(splits[[1]][2]=="quadBM") {
      tmp_file <- paste0("popgrowth_stanmcmc_", file_split2[[1]][1], "_", file_split1[[1]][2], ".RDS")
    }
    file_to_get <- paste0(tmp_dir, "/fits/", tmp_file)
    tmp_mcmc <- readRDS(file_to_get)
    gprint <- ggplot(subset(tmp_mcmc, Parameter == as.character(bad_rhats[i,"parameter"])))+
      geom_line(aes(x=Iteration, y=value, color=as.factor(Chain)))+
      ggtitle(paste(tmp_file, as.character(bad_rhats[i,"parameter"])))
    print(gprint)
  }
  
  if(length(grep("noclimate", as.character(bad_rhats[i,"file"])))==1) {
    splits <- strsplit(as.character(bad_rhats[i,"file"]), "/")
    tmp_dir <- paste0(splits[[1]][1:5], collapse = "/")
    file_split1 <- strsplit(splits[[1]][7],"_")
    file_split2 <- strsplit(file_split1[[1]][4], "[.]")
    if(splits[[1]][2]=="ipm") {
      tmp_file <- paste0("recruitment_stanmcmc_noclimate_", file_split2[[1]][1], "_", file_split1[[1]][2], ".RDS")
    }
    if(splits[[1]][2]=="quadBM") {
      tmp_file <- paste0("popgrowth_stanmcmc_noclimate_", file_split2[[1]][1], "_", file_split1[[1]][2], ".RDS")
    }
    file_to_get <- paste0(tmp_dir, "/fits_noclimate/", tmp_file)
    tmp_mcmc <- readRDS(file_to_get)
    gprint <- ggplot(subset(tmp_mcmc, Parameter == as.character(bad_rhats[i,"parameter"])))+
      geom_line(aes(x=Iteration, y=value, color=as.factor(Chain)))+
      ggtitle(paste(tmp_file, as.character(bad_rhats[i,"parameter"])))
    print(gprint)
  }
}
dev.off()

