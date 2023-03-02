#!/bin/bash
export OMP_NUM_THREADS=4

main_file=$1

for file in `ls $main_file`;do
for inp in `ls $main_file/$file/*.inp`;do
cd $main_file/$file/
echo $main_file/$file/
ls -l *.sh
#sbatch --nodes=1 --ntasks-per-node=4 --job-name=$file 2run.sh
cd ../..
done
done
