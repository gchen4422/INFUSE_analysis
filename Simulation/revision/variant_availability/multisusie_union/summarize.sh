#!/bin/bash
#SBATCH -A pdrineas
#SBATCH -p cpu
#SBATCH --job-name=ms_union_summary
#SBATCH --cpus-per-task=2
#SBATCH --mem=48000
#SBATCH --time=02:00:00
# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -euo pipefail
RUN_DIR=${INFUSE_REVISION}/INFUSE_variant_availability/simulation_native_assv/multisusie_union_test_20261020
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
export TMPDIR="${RUN_DIR}/work/summary_${SLURM_JOB_ID}"
mkdir -p "${TMPDIR}"
Rscript "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/variant_availability/multisusie_union/summarize.R"
