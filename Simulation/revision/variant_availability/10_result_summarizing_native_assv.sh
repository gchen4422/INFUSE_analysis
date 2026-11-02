#!/bin/bash
#SBATCH -A pdrineas
#SBATCH -p cpu
#SBATCH --time=03:00:00
#SBATCH --job-name=native_sum5
#SBATCH --mem=64000
#SBATCH --cpus-per-task=2
#SBATCH --array=0-2

# ---------------------------------------------------------------------------
# Native manuscript summary: ONE ARRAY TASK PER CAUSAL AS-SV LEVEL.
#
# The three levels are independent; each task reads 600 native RDS cells and
# writes only its own results/figures subfolder.
#
#   sbatch 10_result_summarizing_native_assv.sh
# ---------------------------------------------------------------------------
# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
NATIVE_DIR=${INFUSE_REVISION}/INFUSE_variant_availability/simulation_native_assv
mkdir -p "${NATIVE_DIR}/logs_native_sim"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

ASSV_EACH=${SLURM_ARRAY_TASK_ID}

echo "host=$(hostname) task=${SLURM_ARRAY_TASK_ID} native_assv_each=${ASSV_EACH}"
date

Rscript "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/variant_availability/10_result_summarizing_native_assv.R" "${ASSV_EACH}"
status=$?

echo "--- produced for native_ASSV_each_${ASSV_EACH} ---"
REVISION_DIR=${INFUSE_REVISION}
ls -la "${REVISION_DIR}/figures/INFUSE_variant_availability/native_ASSV_each_${ASSV_EACH}" 2>/dev/null | tail -n +2
ls -la "${REVISION_DIR}/res/INFUSE_variant_availability/native_ASSV_each_${ASSV_EACH}" 2>/dev/null | tail -n +2
date
echo "exit ${status}"
exit ${status}
