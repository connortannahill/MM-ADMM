njobs=13
for i in {1..7}
do
  sbatch meshjob$i.sh
done