#!/bin/bash
#SBATCH -A pdrineas
#SBATCH --time=144:00:00
#SBATCH --job-name=simulation
#SBATCH --mem=8000
#SBATCH --array=1-100
#SBATCH --output=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/out/slurm-%A_%a.out
#SBATCH --error=/scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/out/slurm-%A_%a.err


# Get the SLURM task ID (1-based index)
GENE=${SLURM_ARRAY_TASK_ID}

# Define ranges for h2 and causal
H2_VALUES=(1 2)
CAUSAL_VALUES=(1 2 3)
RHO_VALUES=(1 2)



# Loop over h2 and causal combinations for the current gene

for H2 in "${H2_VALUES[@]}"; do
  for CAUSAL in "${CAUSAL_VALUES[@]}"; do
    for RHO in "${RHO_VALUES[@]}"; do
        Rscript --verbose /scratch/negishi/chen4422/hapnest/Simulation_update/simu_region_1mb/shared_50_mesusie_inf/More_annotes/4_shared_50_mesusie_all_annotes_susie_only_ld_outliers_prep.R "${GENE}" "${H2}" "${CAUSAL}" "${RHO}"
    done
  done
done




