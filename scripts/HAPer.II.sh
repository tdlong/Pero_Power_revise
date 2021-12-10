#!/bin/bash
#SBATCH --job-name=simgeno
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=3
#SBATCH --array=1-322

module load R/3.6.2
module load samtools/1.10

file="Phil.parts.HAP.txt"
part=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f1`

Rscript HAPer.II.R "PhilscanHAP" $part

