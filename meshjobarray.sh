#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=def-wan
#SBATCH --mem=50G
#SBATCH --cpus-per-task=32
#SBATCH --array=1-12
#SBATCH --mem=50G
#SBATCH --mail-user=connor.tannahill@uwaterloo.ca
#SBATCH --mail-type=ALL
module load python/3.8.10
module load scipy-stack
python experiments.py < ./SlurmInputFiles/scaletestMon$SLURM_ARRAY_TASK_ID.txt
