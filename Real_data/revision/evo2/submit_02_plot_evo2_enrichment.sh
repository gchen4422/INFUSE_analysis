#!/bin/bash
#SBATCH --job-name=evo2_raw_plot
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00

# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1


Rscript "${INFUSE_ANALYSIS_ROOT}/Real_data/revision/evo2/02_plot_enrichment_either_bp_evo2.R"
Rscript "${INFUSE_ANALYSIS_ROOT}/Real_data/revision/evo2/03_plot_enrichment_either_bmi_evo2.R"
Rscript "${INFUSE_ANALYSIS_ROOT}/Real_data/revision/evo2/04_plot_bmi_asinh_sensitivity.R"
