#!/bin/bash
#SBATCH --job-name=simgeno
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=2
#SBATCH --array=1-338

module load R/3.6.2
module load samtools/1.10

file="Phil.REAL.parts.SNP.txt"
part=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f1`

Rscript SNPer.BT.R "PhilRealSNP" $part

