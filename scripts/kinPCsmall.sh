#!/bin/bash
#SBATCH --job-name=kinPCsmall
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=8

module load R/3.6.2
Rscript kinPCsmall.R

