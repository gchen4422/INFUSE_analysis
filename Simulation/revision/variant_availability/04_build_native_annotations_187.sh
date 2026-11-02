#!/bin/bash
#SBATCH -A pdrineas
#SBATCH -p cpu
#SBATCH --time=01:00:00
#SBATCH --job-name=native_ann187
#SBATCH --mem=8000
#SBATCH --cpus-per-task=1
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
mkdir -p "${NATIVE_DIR}/logs_native_sim" "${NATIVE_DIR}/metadata/annotation_187_qc"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/variant_availability/04_build_native_annotations_187.R" "${SLURM_ARRAY_TASK_ID}"
status=$?
exit ${status}
