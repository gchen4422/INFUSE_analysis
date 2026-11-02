#!/bin/bash
# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -euo pipefail

# Submit four scenario arrays, then summarize and plot only after every array
# completes successfully. Existing validated cell results are skipped safely.

SCENARIOS=(baseline different_annotation different_weight combined)
JOB_IDS=()

for SCENARIO in "${SCENARIOS[@]}"; do
  JOB_ID=$(sbatch --parsable \
    --job-name="paintor_${SCENARIO}" \
    "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/annotation_heterogeneity/06_run_paintor_benchmarks_cell.sh" "${SCENARIO}")
  JOB_IDS+=("${JOB_ID}")
  echo "Submitted ${SCENARIO}: ${JOB_ID}"
done

DEPENDENCY=$(IFS=:; echo "${JOB_IDS[*]}")
SUMMARY_JOB=$(sbatch --parsable \
  --dependency="afterok:${DEPENDENCY}" \
  "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/annotation_heterogeneity/07_summarize_paintor_benchmarks.sh")
PLOT_JOB=$(sbatch --parsable \
  --dependency="afterok:${SUMMARY_JOB}" \
  "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/annotation_heterogeneity/08_plot_combined_with_paintor.sh")

echo "Submitted summary: ${SUMMARY_JOB}"
echo "Submitted combined plot: ${PLOT_JOB}"
