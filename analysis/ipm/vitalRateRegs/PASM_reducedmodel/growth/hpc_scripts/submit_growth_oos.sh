#!/bin/bash
#SBATCH --array=1-312
#SBATCH --job-name=growth_oos
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4-00:00:00

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH -$SLURM_ARRAY_TASK_ID $HOME/growth_oos/growth_oos_launcher.R