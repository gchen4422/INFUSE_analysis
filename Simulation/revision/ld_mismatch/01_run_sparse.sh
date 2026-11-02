#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --array=1-500
# Set these paths in the environment before running.
: "${INFUSE_ANALYSIS_ROOT:?Set INFUSE_ANALYSIS_ROOT to the analysis repository}"
: "${INFUSE_PROJECT:?Set INFUSE_PROJECT to the baseline simulation data}"
: "${INFUSE_REVISION:?Set INFUSE_REVISION to a separate output directory}"
export INFUSE_ANALYSIS_ROOT INFUSE_PROJECT INFUSE_REVISION
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
if [ -n "${INFUSE_R_LIBRARY:-}" ]; then export R_LIBS_USER="${INFUSE_R_LIBRARY}${R_LIBS_USER:+:${R_LIBS_USER}}"; fi
set -euo pipefail
cd "${INFUSE_PROJECT}"
module_dir="${INFUSE_REVISION}/LD_mismatch_sensitivity"
mkdir -p "${module_dir}/logs" "${module_dir}/runtime"
export TMPDIR="${module_dir}/runtime"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
scenarios=(matched sparse_attenuation_025 sparse_attenuation_050 sparse_negative_050 sparse_negative_100)
task=${SLURM_ARRAY_TASK_ID:?}
if (( task < 1 || task > 500 )); then exit 1; fi
scenario_index=$(( (task - 1) / 100 ))
locus=$(( (task - 1) % 100 + 1 ))
Rscript --vanilla "${INFUSE_ANALYSIS_ROOT}/Simulation/revision/ld_mismatch/01_run_sparse.R" "${locus}" "${scenarios[$scenario_index]}"
