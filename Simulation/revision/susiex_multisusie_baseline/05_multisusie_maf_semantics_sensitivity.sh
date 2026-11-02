#!/bin/bash
#SBATCH -A pdrineas
#SBATCH -p cpu
#SBATCH --time=02:00:00
#SBATCH --job-name=base_maf_check
#SBATCH --mem=8000
#SBATCH --cpus-per-task=2

# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -u
MODULE_DIR=${INFUSE_REVISION}/baseline_susiex_multisusie_benchmark
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX="${MODULE_DIR}/.pycache"
export RETICULATE_PYTHON=${INFUSE_PYTHON_BIN}
Rscript "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/susiex_multisusie_baseline/05_multisusie_maf_semantics_sensitivity.R"

