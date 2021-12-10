#!/bin/bash
#SBATCH --job-name=simgeno
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=2
 
module load R/3.6.2
module load samtools/1.10

## consolidate
cat simgeno/HAP_1_Chr23.txt simgeno/HAP_2_Chr23.txt simgeno/HAP_3_Chr23.txt simgeno/HAP_4_Chr23.txt simgeno/HAP_5_Chr23.txt\
    simgeno/HAP_6_Chr23.txt simgeno/HAP_7_Chr23.txt simgeno/HAP_8_Chr23.txt simgeno/HAP_9_Chr23.txt simgeno/HAP_10_Chr23.txt\
    simgeno/HAP_11_Chr23.txt simgeno/HAP_12_Chr23.txt simgeno/HAP_13_Chr23.txt simgeno/HAP_14_Chr23.txt simgeno/HAP_15_Chr23.txt\
    simgeno/HAP_16_Chr23.txt simgeno/HAP_17_Chr23.txt simgeno/HAP_18_Chr23.txt simgeno/HAP_19_Chr23.txt simgeno/HAP_20_Chr23.txt\
    | tr [:blank:] \\t | gzip -c > simgeno/HAP_Chr23.txt.gz 
paste simgeno/SNP_1_Chr23.txt simgeno/SNP_2_Chr23.txt simgeno/SNP_3_Chr23.txt simgeno/SNP_4_Chr23.txt simgeno/SNP_5_Chr23.txt\
      simgeno/SNP_6_Chr23.txt simgeno/SNP_7_Chr23.txt simgeno/SNP_8_Chr23.txt simgeno/SNP_9_Chr23.txt simgeno/SNP_10_Chr23.txt\
      simgeno/SNP_11_Chr23.txt simgeno/SNP_12_Chr23.txt simgeno/SNP_13_Chr23.txt simgeno/SNP_14_Chr23.txt simgeno/SNP_15_Chr23.txt\
      simgeno/SNP_16_Chr23.txt simgeno/SNP_17_Chr23.txt simgeno/SNP_18_Chr23.txt simgeno/SNP_19_Chr23.txt simgeno/SNP_20_Chr23.txt\
    | tr [:blank:] \\t | gzip -c > simgeno/SNP_Chr23.txt.gz 
rm simgeno/HAP_*_Chr23.txt
rm simgeno/SNP_*_Chr23.txt

