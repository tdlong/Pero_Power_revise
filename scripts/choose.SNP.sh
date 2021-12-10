#!/bin/bash
#SBATCH --job-name=pickSNPs
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=4  

module load R/3.2.1
Rscript choose.SNP.R

