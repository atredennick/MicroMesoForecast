#!/bin/bash
sbatch ./ipm_growth_validation/submit_growth_loyo.sh
sbatch ./ipm_surv_validation/submit_survival_BOGR_loyo.sh
sbatch ./ipm_surv_validation/submit_survival_HECO_loyo.sh
sbatch ./ipm_surv_validation/submit_survival_PASM_loyo.sh
sbatch ./ipm_surv_validation/submit_survival_POSE_loyo.sh
sbatch ./ipm_rec_validation/submit_recruitment_loyo.sh