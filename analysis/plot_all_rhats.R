##  Script to plot all r-hats for vital rate model parameters

rm(list=ls(all.names = TRUE))

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
