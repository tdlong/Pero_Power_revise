#!/bin/bash
#SBATCH --job-name=fakepower
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=1
#SBATCH --array=1-320
 
module load R/3.6.2
module load samtools/1.10

file="Chr23.hunks.HAP.txt"
Nind=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f1`
part=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f2`
Rscript fake.power.HAP.R "simgeno" $Nind $part

