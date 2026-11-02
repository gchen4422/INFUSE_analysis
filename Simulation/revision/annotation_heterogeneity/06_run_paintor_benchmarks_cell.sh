#!/bin/bash
#SBATCH --job-name=paintor_prior
#SBATCH --array=1-600
#SBATCH --mem=6G
#SBATCH --time=48:00:00
# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -euo pipefail

# Usage:
#   sbatch 06_run_paintor_benchmarks_cell.sh <scenario> [task_id]
#
# One 600-cell array covers 100 loci x 2 h2 settings x 3 causal-count
# settings for one scenario. Each cell runs Paintor, Paintor+ (the same 12
# columns used by INFUSE), and Paintor+H (the same 189 used by INFUSE-H).

LOG_DIR="${INFUSE_REVISION}/ancestry_specific_eaf_prior_simulation/logs/paintor"
SCENARIO="${1:?scenario is required}"
TASK_ID="${2:-${SLURM_ARRAY_TASK_ID:-1}}"

mkdir -p "${LOG_DIR}"

if (( TASK_ID < 1 || TASK_ID > 600 )); then
  echo "TASK_ID must be between 1 and 600." >&2
  exit 2
fi

ZERO_BASED=$((TASK_ID - 1))
LOCUS=$((ZERO_BASED % 100 + 1))
COMBINATION=$((ZERO_BASED / 100))
H2_INDEX=$((COMBINATION % 2 + 1))
CAUSAL_INDEX=$((COMBINATION / 2 + 1))

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

Rscript --vanilla \
  "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/annotation_heterogeneity/06_run_paintor_benchmarks_cell.R" \
  "${SCENARIO}" "${LOCUS}" "${H2_INDEX}" "${CAUSAL_INDEX}"
