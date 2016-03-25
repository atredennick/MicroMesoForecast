#!/bin/bash
#SBATCH --array=1-312
#SBATCH --job-name=surv_oos
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --time=4-00:00:00
#SBATCH --constraint=gonium

# Send mail
#SBATCH --mail-type=ALL
#SBATCH --mail-user=atredenn@gmail.com

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH -$SLURM_ARRAY_TASK_ID $HOME/survival_oos/survival_oos_launcher.R