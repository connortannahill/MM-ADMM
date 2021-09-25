#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --account=def-wan
#SBATCH --cpus-per-task=32
python experiments.py < scaletest.txt
