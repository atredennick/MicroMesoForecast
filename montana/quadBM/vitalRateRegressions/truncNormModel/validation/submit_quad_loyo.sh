#!/bin/bash
#SBATCH --job-name=R_quadLOYO_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4-00:00:00
#SBATCH --constraint=gonium

# Send mail
#SBATCH --mail-type=ALL
#SBATCH --mail-user=atredenn@gmail.com

#SBATCH --array=1933-1945

. /rc/tools/utils/dkinit
use JAGS
reuse -q R

R CMD BATCH $HOME/quadmodel/growthAllSppHierarchy_MCMC_LOYO_$SLURM_ARRAY_TASK_ID.R