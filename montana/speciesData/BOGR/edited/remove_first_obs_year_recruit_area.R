##  This script subsets out recruit area for quadrats that no longer
##    have recruit areas that can be fit because we removed the years
##    before that event. We took out quad-years for BOGR that looked
##    suspicious in terms of mapping. Now that we took those out, we
##    need to also update the recruitArea.csv file so that we do not
##    have recruit areas for what is now the first year of our observations.

##  Author:       Andrew Tredennick
##  Email:        atredenn@gmail.com
##  Last update:  4-21-2015

##  Note that this script requires a connection the Adler Lab USU server
##    for access to files that I don't want to make multiple copies of.

# Clear the workspace 
rm(list=ls(all=TRUE))

# Load libraries
library(reshape2)
library(plyr)

####
####  Read in necessary data files ---------------------------------------------
####
bad_bogrs <- read.csv("/Volumes/adlerlab/group/montanachart/suspect_BOGR_quads.csv")
current_recruit_area <- read.csv("recArea_quads_not_removed.csv")


####
####  Find instances where initial year(s) were removed ------------------------
####
# Find minimum year for each quad
bad_bogrs_beginyear <- ddply(bad_bogrs, .(quadrat), summarise,
                             first_year = min(year))
current_beginyear <- ddply(current_recruit_area, .(quad), summarise,
                           first_year = min(year, na.rm = TRUE))

# Look at coincident quads and remove recruit minyear rows 
#   IF minyear from current != minyear from bad_bogrs
badquads <- which(current_beginyear$quad %in% bad_bogrs_beginyear$quadrat)
current_badquads_only <- current_beginyear[badquads,]
bad_bogrs_beginyear$first_year <- bad_bogrs_beginyear$first_year - 1900
first_years_to_remove <- which(current_badquads_only$first_year > bad_bogrs_beginyear$first_year)
current_badquads_only$remove_flag <- "keep"
current_badquads_only$remove_flag[first_years_to_remove] <- "remove"

# Merge the current bad quads back in with full inventory
current_inventory_tmp <- merge(current_beginyear, current_badquads_only, all=TRUE)
current_inventory <- merge(current_inventory_tmp, current_recruit_area,
                           by.x=c("quad", "first_year"), by.y=c("quad", "year"),
                           all=TRUE)
rows_to_remove <- which(current_inventory$remove_flag == "remove")
final_inventory <- current_inventory[-rows_to_remove,]


####
####  Write to file ------------------------------------------------------------
####
col_remove <- which(colnames(final_inventory) %in% c("remove_flag", "X"))
final_inventory <- final_inventory[,-col_remove]
cols_to_name <- which(colnames(current_recruit_area) != "X")
colnames(final_inventory) <- colnames(current_recruit_area)[cols_to_name]
write.csv(final_inventory, "recArea.csv")

