#!/bin/bash
#SBATCH --job-name=SNPer
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=2
#SBATCH --array=11-322

module load R/3.6.2
module load samtools/1.10

file="Phil.parts.SNP.txt"
part=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f1`

Rscript SNPer.sim4500.R "PhilscanSNP" $part

