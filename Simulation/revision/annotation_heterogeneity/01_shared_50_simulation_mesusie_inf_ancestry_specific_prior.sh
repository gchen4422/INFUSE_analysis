#!/bin/bash
#SBATCH --job-name=infuse_asp
#SBATCH --array=1-100
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -euo pipefail

# Usage:
#   bash 01_shared_50_simulation_mesusie_inf_ancestry_specific_prior.sh \
#     <scenario> <base_d> [LD_BLOCK]
#
# Scenarios:
#   baseline              same annotations, weights (2, 2)
#   different_annotation  ancestry-specific aligned EAF, weights (2, 2)
#   different_weight      same annotations, weights (2, 1)
#   combined              ancestry-specific aligned EAF, weights (2, 1)
#
# base_d is 10 or 187. EAF_EUR and EAF_AFR are appended, so the annotation
# matrix supplied downstream contains 12 or 189 columns.
# When submitted with sbatch, each array task handles one LD block. When run
# directly with bash, LD_BLOCK can be supplied explicitly; if it is omitted,
# the R script processes loci 1 through 100 sequentially.
#
# Examples:
#   sbatch 01_shared_50_simulation_mesusie_inf_ancestry_specific_prior.sh \
#     different_annotation 10
#   bash 01_shared_50_simulation_mesusie_inf_ancestry_specific_prior.sh \
#     different_annotation 10 1

SCENARIO="${1:-different_annotation}"
BASE_D="${2:-10}"
LD_BLOCK="${3:-${SLURM_ARRAY_TASK_ID:-}}"

Rscript "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/annotation_heterogeneity/01_shared_50_simulation_mesusie_inf_ancestry_specific_prior.R" \
  "${SCENARIO}" "${BASE_D}" "${LD_BLOCK}"
