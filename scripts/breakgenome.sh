#!/bin/bash
#SBATCH --job-name=breakgenome
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=2

rm simgeno/HAP.all.sort.txt
rm simgeno/SNP.all.sort.txt
KKK=(1 2 3 4 5 6 7 8a 8b 9 10 11 12 13 14 15 1621 17 18 20 22 23)
for i in "${KKK[@]}"
do
	zcat PhilsampsII/HAPs.Chr${i}.censor.rf.txt.gz >> PhilsampsII/HAP.all.sort.txt
	zcat PhilsampsII/SNPs.Chr${i}.censor.txt.gz >> PhilsampsII/SNP.all.sort.txt
done
# 348*5000 = 1740000
cat PhilsampsII/HAP.all.sort.txt | split - -l 1740000 --filter='gzip > $FILE.gz' simgenoHunksII/HAPs.part.
cat PhilsampsII/SNP.all.sort.txt | split - -l 50000 --filter='gzip > $FILE.gz' simgenoHunksII/SNPs.part.
## Background is every 400th SNP (starting at 200)
## So I want to make sure the SNPs used to contruct the relatedness matrix are "out of phase" with background SNP
cat PhilsampsII/SNP.all.sort.txt | awk '{if(NR % 631 == 137) print}' | awk '{for(i=3;i<NF;i++){printf("%d\t",$i+0.5)};{i=NF;printf("%d\n",$i+0.5)}}' >hard.calls.txt

ls simgenoHunksII/SNPs.part* | cut -f2 -d"/" | cut -f3 -d"." > Phil.parts.SNP.txt   #322
ls simgenoHunksII/HAPs.part* | cut -f2 -d"/" | cut -f3 -d"." > Phil.parts.HAP.txt   #322

