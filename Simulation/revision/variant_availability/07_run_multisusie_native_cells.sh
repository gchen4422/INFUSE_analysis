#!/bin/bash
#SBATCH -A pdrineas
#SBATCH -p cpu
#SBATCH --time=04:00:00
#SBATCH --job-name=native_multisusie
#SBATCH --mem=8000
#SBATCH --cpus-per-task=2
#SBATCH --array=1-900

# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -u
NATIVE_DIR=${INFUSE_REVISION}/INFUSE_variant_availability/simulation_native_assv
mkdir -p "${NATIVE_DIR}/logs_native_sim" "${NATIVE_DIR}/results_native_sim/multisusie_intersect"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

TASK_OFFSET=${TASK_OFFSET:-0}
global_task_id=$((SLURM_ARRAY_TASK_ID + TASK_OFFSET))
x=$((global_task_id - 1))
locus=$((x % 100 + 1)); x=$((x / 100))
base_index=$((x % 3 + 1)); x=$((x / 3))
h2_index=$((x % 2 + 1)); x=$((x / 2))
n_assv_each=$((x % 3))

echo "host=$(hostname) local_task=${SLURM_ARRAY_TASK_ID} global_task=${global_task_id} locus=${locus} base_index=${base_index} h2_index=${h2_index} n_assv_each=${n_assv_each}"
date
Rscript "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/variant_availability/07_run_multisusie_native_cell.R" "${locus}" \
  "${base_index}" "${h2_index}" "${n_assv_each}"
status=$?
date
echo "exit ${status}"
exit ${status}
