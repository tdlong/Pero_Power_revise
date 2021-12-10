#!/bin/bash
#SBATCH --job-name=SNPerPro
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=1
#SBATCH --array=11-322

file="Phil.parts.SNP.txt"
part=`head -n $SLURM_ARRAY_TASK_ID $file | tail -n 1 | cut -f1`

python findtop.new.LOD7.5.py "PhilscanSNP" $part

