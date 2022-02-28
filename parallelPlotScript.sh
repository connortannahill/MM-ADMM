#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-wan
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=connor.tannahill@uwaterloo.ca
#SBATCH --mail-type=ALL

module load python/3.8.10
module load scipy-stack

monNums=( 1 2 3)

for monNum in "${monNums[@]}"; do
echo $monNum

# inputStr="plot_energy_decrease()\nMonitor"$num"\nTrue\nTrue\n"
# echo $inputStr
python experiments.py <<STDIN -o other --options
create_parallel_plot()
3DMonitor$monNum
False
exit()
STDIN

done
