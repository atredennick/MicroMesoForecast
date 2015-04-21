Vital rates based on a subset of BOGR quadrats that excluded the "fishy" ones. 
See H:\montanachart\lifetables\polys\Species\BOGR-edited

ATT: recArea_quads_not_removed.csv is the recruitment area data with first non-fishy BOGR quad years still present. We need to remove these since they are now the beginning year with genets age=1. The script remove_first_obs_year_recruit_area.R does this and creates a new file called recArea.csv without those age=1 recruit areas based on the quad inventory with bad BOGRs removed. That is the file to use in all subsequent analysis.