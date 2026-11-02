#!/bin/bash
# Submit extraction for d=10 and d=187, then analyze and plot after success.
# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -euo pipefail

MODULE="${INFUSE_REVISION}/configuration_classification_w2"
LOG_DIR="${MODULE}/logs"
mkdir -p "${LOG_DIR}"

extract_job=$(sbatch --parsable \
  --job-name=config_extract \
  --array=1-2 \
  --cpus-per-task=1 --mem=8G --time=4:00:00 \
  --output="${LOG_DIR}/extract_%A_%a.out" \
  --error="${LOG_DIR}/extract_%A_%a.err" \
  --wrap="export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1; label=\$( [ \${SLURM_ARRAY_TASK_ID} -eq 1 ] && echo d10 || echo d187 ); Rscript '${INFUSE_ANALYSIS_ROOT}/Simulation/revision/configuration_classification/01_extract_configuration_predictions.R' \${label}")

analysis_job=$(sbatch --parsable \
  --job-name=config_analyze \
  --dependency="afterok:${extract_job}_1:${extract_job}_2" \
  --cpus-per-task=1 --mem=8G --time=4:00:00 \
  --output="${LOG_DIR}/analyze_%j.out" \
  --error="${LOG_DIR}/analyze_%j.err" \
  --wrap="export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1; Rscript '${INFUSE_ANALYSIS_ROOT}/Simulation/revision/configuration_classification/02_analyze_configuration_classification.R'; Rscript '${INFUSE_ANALYSIS_ROOT}/Simulation/revision/configuration_classification/03_plot_configuration_classification.R'; Rscript '${INFUSE_ANALYSIS_ROOT}/Simulation/revision/configuration_classification/04_analyze_credible_set_configuration.R'; Rscript '${INFUSE_ANALYSIS_ROOT}/Simulation/revision/configuration_classification/05_plot_credible_set_configuration.R'")

printf 'Extraction array: %s\nAnalysis/plot job: %s\n' "${extract_job}" "${analysis_job}"
