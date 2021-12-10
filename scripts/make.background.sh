#!/bin/bash
#SBATCH --job-name=background
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --array=1-23
#SBATCH --cpus-per-task=1  

# skip X...
KKK=(1 2 3 4 5 6 7 8a 8b 9 10 11 12 13 14 15 1621 17 18 19 20 22 23 X)
mychr=${KKK[$SLURM_ARRAY_TASK_ID-1]}

dir="PhilsampsII"

# just cumulative sum of every 400th SNP
zcat $dir/SNPs.Chr${mychr}.censor.txt.gz | cut -d " " -f 3- |\
    awk '{if(NR%400==200){print $0}}' |\
    awk '{ for (i=1; i<=NF; i++) total[i] += $i; }; END { for (i=1; i<=NF; i++) {printf "%s ",total[i]};print "" }' >\
    $dir/SNPs.Chr${mychr}.BG.txt            

