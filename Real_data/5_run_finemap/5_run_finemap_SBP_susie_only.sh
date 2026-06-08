#!/bin/bash
#SBATCH --array=1-46    # Define the SLURM array range
#SBATCH --time=24:00:00
#SBATCH -A pdrineas
#SBATCH --mem=16000
#SBATCH --job-name=run_finemap_susie
#SBATCH --output=logs/run_finemap_%A_%a.out
#SBATCH --error=logs/run_finemap_%A_%a.err

# SLURM provides the task ID for the array job in the $SLURM_ARRAY_TASK_ID variable.
sbatch_array_ID=${SLURM_ARRAY_TASK_ID}

# Run the R script with the array ID as an argument
Rscript 5_run_finemap_SBP_susie_only.R ${sbatch_array_ID}
