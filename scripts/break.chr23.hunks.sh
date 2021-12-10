#!/bin/bash
#SBATCH --job-name=break
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=2
zcat simgeno/HAP_Chr23.txt.gz | tail -n +2 | sort -k2,2n -k3,3n > simgeno/HAP_Chr23.sort.txt
zcat simgeno/SNP_Chr23.txt.gz | tail -n +2 | split - -l 10000 --filter='gzip > $FILE.gz' simgenoHunks/SNPs.part.
cat simgeno/HAP_Chr23.sort.txt | split - -l 1780000 --filter='gzip > $FILE.gz' simgenoHunks/HAPs.part.

