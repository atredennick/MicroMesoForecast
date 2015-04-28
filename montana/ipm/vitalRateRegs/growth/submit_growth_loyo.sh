#!/bin/bash
#SBATCH --array=32-45
#SBATCH --job-name=R_growthLOYO_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4-00:00:00
#SBATCH --constraint=gonium

# Send mail
#SBATCH --mail-type=ALL
#SBATCH --mail-user=atredenn@gmail.com

. /rc/tools/utils/dkinit
use JAGS
reuse -q R

R CMD BATCH -$SLURM_ARRAY_TASK_ID $HOME/montana_vitalrates_ipm/growthAllSpp_MCMC_LOYO.R