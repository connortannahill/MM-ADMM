#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-wan
#SBATCH --mem=50G
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=connor.tannahill@uwaterloo.ca
#SBATCH --mail-type=ALL
module load python/3.8.10
module load scipy-stack
python experiments.py < scaletestMon8.txt
