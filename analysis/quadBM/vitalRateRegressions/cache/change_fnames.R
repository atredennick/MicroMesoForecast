# Quick and dirty script to rename QBM validation fits to get years lined up with new version
# Just takes the name and subtracts 1 year from the leavout year id

setwd("/Users/atredenn/Repos/MicroMesoForecast/analysis/quadBM/vitalRateRegressions/truncNormModel/validation/fits/")
all_files <- list.files()
for(ifile in 1:length(all_files)){
  tmpchar <- all_files[ifile]
  oldyr <- strsplit(strsplit(tmpchar,"out")[[1]][2], "[.]")[[1]][1]
  newyr <- as.character(as.numeric(oldyr)-1)
  new_fname <- paste0(strsplit(tmpchar,"out")[[1]][1],"out",newyr,".RDS")
  file.rename(tmpchar, new_fname)
}