#!/bin/bash
#SBATCH -A pdrineas
#SBATCH -p cpu
#SBATCH --time=01:00:00
#SBATCH --job-name=w2_pip_cal
#SBATCH --mem=32000
#SBATCH --cpus-per-task=2
#SBATCH --array=1-2

# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
ANALYSIS=${INFUSE_REVISION}/pip_calibration_w2_d10_d187
mkdir -p "${ANALYSIS}/logs"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

case "${SLURM_ARRAY_TASK_ID}" in
  1) setting=d10 ;;
  2) setting=d187 ;;
  *) echo "invalid array task ${SLURM_ARRAY_TASK_ID}" >&2; exit 2 ;;
esac

echo "host=$(hostname) setting=${setting}"
date
Rscript "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/pip_calibration/01_generate_calibration.R" "${setting}"
status=$?
date
echo "exit status ${status}"
exit ${status}
