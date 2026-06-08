#!/bin/bash
#SBATCH --job-name=evo2_dbp
#SBATCH --account=ppaschou
#SBATCH --partition=a100-40gb
#SBATCH --gres=gpu:1
#SBATCH --mem=50G
#SBATCH --qos=standby
#SBATCH --time=04:00:00
#SBATCH --array=1-38                   # The magic line: Jobs 1-38
#SBATCH --output=logs/evo2_%A_%a.out     # %A is the master job ID, %a is the array ID
#SBATCH --error=logs/evo2_%A_%a.err      # Standard error log

# 1. Navigate to your working directory
cd /scratch/gilbreth/chen4422/hapnest/to_depot_locus_to_evo2/dbp_summ_stat/

# 2. Load the environment
module load anaconda
source activate /depot/gao824/chen4422/conda_env_list/evo2_env

# 3. Dynamically set the filenames using the SLURM_ARRAY_TASK_ID variable
INPUT_FILE="eur_loci_${SLURM_ARRAY_TASK_ID}_summ"
OUTPUT_FILE="loci_${SLURM_ARRAY_TASK_ID}_summ_scored.csv"

echo "--------------------------------------------------------"
echo "Starting task ${SLURM_ARRAY_TASK_ID} on node ${HOSTNAME}"
echo "Input: ${INPUT_FILE} | Output: ${OUTPUT_FILE}"
echo "--------------------------------------------------------"

# 4. Run the python script
python /scratch/gilbreth/chen4422/hapnest/to_depot_locus_to_evo2/get_evo2_score/score_evo2.py --input ${INPUT_FILE} --output ${OUTPUT_FILE}
