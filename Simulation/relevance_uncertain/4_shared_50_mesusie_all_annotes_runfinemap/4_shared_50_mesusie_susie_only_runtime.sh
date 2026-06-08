#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=4:00:00
#SBATCH --job-name=simulation
#SBATCH --mem=16000
#SBATCH --array=1-100
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/out/slurm-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/out/slurm-%A_%a.err
#SBATCH --partition=cpu
#SBATCH --qos=standby


set -euo pipefail

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Get the SLURM task ID (1-based index)
GENE=${SLURM_ARRAY_TASK_ID}

# Define ranges for h2 and causal
H2_VALUES=(1 2)
CAUSAL_VALUES=(1 2 3)

# Loop over h2 and causal combinations for the current gene
for H2 in "${H2_VALUES[@]}"; do
  for CAUSAL in "${CAUSAL_VALUES[@]}"; do
    # Run the R script for each combination
    Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/4_shared_50_mesusie_susie_only_runtime.R "${GENE}" "${H2}" "${CAUSAL}"
  done
done







