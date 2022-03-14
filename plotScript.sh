#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --account=def-wan
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=connor.tannahill@uwaterloo.ca
#SBATCH --mail-type=ALL

module load python/3.8.10
module load scipy-stack

# nums=(10, 20, 40, 80, 160, 320)
monNums=( 1 2 3)
# monNums=( 3)
nums=( 10 20 40 80 160 320)
# echo "${nums[@]}"

for num in "${nums[@]}"; do
for monNum in "${monNums[@]}"; do

echo $num
echo $monNum
# inputStr="plot_energy_decrease()\nMonitor"$num"\nTrue\nTrue\n"
# echo $inputStr
python experiments.py <<STDIN -o other --options
plot_energy_decrease()
3DMonitor$monNum$num
True
True
exit()
STDIN

done
done
