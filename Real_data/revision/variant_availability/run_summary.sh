#!/bin/bash
# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -euo pipefail
bash "${INFUSE_ANALYSIS_ROOT}/Real_data/revision/variant_availability/run_job.sh" 08_summarize_results.R
bash "${INFUSE_ANALYSIS_ROOT}/Real_data/revision/variant_availability/run_job.sh" 09_matched_size_enrichment.R
bash "${INFUSE_ANALYSIS_ROOT}/Real_data/revision/variant_availability/run_job.sh" 10_pooled_configuration_enrichment.R
bash "${INFUSE_ANALYSIS_ROOT}/Real_data/revision/variant_availability/run_job.sh" 11_plot_enrichment.R
bash "${INFUSE_ANALYSIS_ROOT}/Real_data/revision/variant_availability/run_job.sh" 12_finalize_matched_summary.R
bash "${INFUSE_ANALYSIS_ROOT}/Real_data/revision/variant_availability/run_job.sh" 13_cs_summary_panels.R
bash "${INFUSE_ANALYSIS_ROOT}/Real_data/revision/variant_availability/run_job.sh" 14_combine_real_data_panels.R
