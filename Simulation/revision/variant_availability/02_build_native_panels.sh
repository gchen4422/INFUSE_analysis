#!/bin/bash
#SBATCH -A pdrineas
#SBATCH -p cpu
#SBATCH --time=04:00:00
#SBATCH --job-name=native_assv_panel
#SBATCH --mem=16000
#SBATCH --cpus-per-task=2
#SBATCH --array=1-100

# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -u
NATIVE_DIR=${INFUSE_REVISION}/INFUSE_variant_availability/simulation_native_assv
mkdir -p "${NATIVE_DIR}/logs_source_preserved" "${NATIVE_DIR}/.work_source_preserved"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export FORCE=${FORCE:-0}

echo "host=$(hostname) locus=${SLURM_ARRAY_TASK_ID} force=${FORCE}"
date
Rscript "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/variant_availability/02_build_native_locus.R" "${SLURM_ARRAY_TASK_ID}"
status=$?
date
echo "exit ${status}"
exit ${status}
