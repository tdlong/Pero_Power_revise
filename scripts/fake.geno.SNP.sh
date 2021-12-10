#!/bin/bash
#SBATCH --job-name=fakepower
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=2
#SBATCH --array=1-288
 
module load R/3.6.2
module load samtools/1.10

file="Chr23.hunks.SNP.txt"
Nind=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f1`
part=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f2`
Rscript fake.power.SNP.R "simgeno" $Nind $part

