#!/bin/bash
#SBATCH --array=1-13
#SBATCH --job-name=R_quadLOYO_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4-00:00:00
#SBATCH --constraint=gonium

# Send mail
#SBATCH --mail-type=ALL
#SBATCH --mail-user=atredenn@gmail.com

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH -$SLURM_ARRAY_TASK_ID $HOME/pop_validation/popgrowth_STAN_validation.R