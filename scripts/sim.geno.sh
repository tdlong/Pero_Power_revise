#!/bin/bash
#SBATCH --job-name=simgeno
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=4
#SBATCH --array=1-20
 
module load R/3.6.2
module load samtools/1.10

Rscript sim.geno.R "Chr23" 200 $SLURM_ARRAY_TASK_ID simgeno

