#!/bin/bash
#SBATCH --array=1-199
#SBATCH --time=4:00:00
#SBATCH -A pdrineas
#SBATCH --mem=32000
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu
#SBATCH --qos=standby
#SBATCH --job-name=run_finemap_susie
#SBATCH --output=logs/run_finemap_%A_%a.out
#SBATCH --error=logs/run_finemap_%A_%a.err

set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

sbatch_array_ID=${SLURM_ARRAY_TASK_ID}

Rscript --vanilla 5_run_finemap_BMI_susie_only_evo2.R "${sbatch_array_ID}"
