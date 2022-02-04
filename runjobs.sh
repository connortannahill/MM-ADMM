njobs=13
for i in {1..13}
do
  sbatch meshjob$i.sh
done