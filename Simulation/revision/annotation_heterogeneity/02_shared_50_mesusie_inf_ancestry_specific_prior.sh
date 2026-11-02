#!/bin/bash
#SBATCH --job-name=infuse_eval
#SBATCH --array=1-600%100
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=16:00:00
# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -euo pipefail

# Usage:
#   sbatch 02_shared_50_mesusie_inf_ancestry_specific_prior.sh \
#     <scenario> <base_d> [result_tag]
#
# Scenarios: baseline, different_annotation, different_weight, combined
# base_d: 10 or 187 (the fitted annotation matrix has 12 or 189 columns)
#
# The 600 array tasks span 100 loci x 2 h2 settings x 3 causal-count
# settings. A fourth argument can provide a single task index for direct testing.

SCENARIO="${1:-different_annotation}"
BASE_D="${2:-10}"
RESULT_TAG="${3:-}"
TASK_ID="${4:-${SLURM_ARRAY_TASK_ID:-1}}"

if (( TASK_ID < 1 || TASK_ID > 600 )); then
  echo "TASK_ID must be between 1 and 600." >&2
  exit 2
fi

ZERO_BASED=$((TASK_ID - 1))
LD_BLOCK=$((ZERO_BASED % 100 + 1))
COMBINATION=$((ZERO_BASED / 100))
H2_NUM=$((COMBINATION % 2 + 1))
NUM_CAUSAL=$((COMBINATION / 2 + 1))

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

Rscript --vanilla \
  "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/annotation_heterogeneity/02_shared_50_mesusie_inf_ancestry_specific_prior.R" \
  "${SCENARIO}" "${BASE_D}" "${LD_BLOCK}" "${H2_NUM}" "${NUM_CAUSAL}" \
  "${RESULT_TAG}"
