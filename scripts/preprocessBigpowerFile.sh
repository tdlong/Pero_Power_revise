#!/bin/bash
#SBATCH --job-name=procbigpow
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=6

module load R/3.6.2
Rscript preprocessBigpowerFile.R

