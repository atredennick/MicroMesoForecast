#!/bin/bash
#SBATCH --array=1-52
#SBATCH --job-name=R_survLOYO_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4-00:00:00

. /rc/tools/utils/dkinit
reuse -q R

R CMD BATCH -$SLURM_ARRAY_TASK_ID $HOME/ipm_surv_validation/survivalAllSpp_validation.R