#!/bin/bash
# Invoke with sbatch and an explicit --chdir, --output and --error in this project.
# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -euo pipefail
ROOT=${AOU_UNION_ROOT:-${INFUSE_REVISION}/INFUSE_variant_availability/real_data_AoU_union}
export AOU_UNION_ROOT="$ROOT"
export TMPDIR="$ROOT/tmp/${SLURM_JOB_ID:-interactive}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMPDIR"
export TMP="$TMPDIR" TEMP="$TMPDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
ulimit -c 0
cd "$ROOT"
case "$1" in
  *.py) "${INFUSE_PYTHON_BIN:-python3}" "${INFUSE_ANALYSIS_ROOT}/Real_data/revision/variant_availability/$1" "${@:2}" ;;
  *) Rscript --vanilla "${INFUSE_ANALYSIS_ROOT}/Real_data/revision/variant_availability/$1" "${@:2}" ;;
esac
