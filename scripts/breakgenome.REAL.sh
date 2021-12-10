#!/bin/bash
#SBATCH --job-name=breakgenome
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --cpus-per-task=2

rm simgeno/HAP.all.sort.txt
rm simgeno/SNP.all.sort.txt
KKK=(1 2 3 4 5 6 7 8a 8b 9 10 11 12 13 14 15 1621 17 18 19 20 22 23)
for i in "${KKK[@]}"
do
	zcat PhilsampsII/HAPs.Chr${i}.censor.rf.txt.gz >> simgenoHunksIII/HAP.all.sort.txt
	zcat PhilsampsII/SNPs.Chr${i}.censor.txt.gz >> simgenoHunksIII/SNP.all.sort.txt
done
# 348*5000 = 1740000
cat simgenoHunksIII/HAP.all.sort.txt | split - -l 1740000 --filter='gzip > $FILE.gz' simgenoHunksIII/HAPs.part.
cat simgenoHunksIII/SNP.all.sort.txt | split - -l 50000 --filter='gzip > $FILE.gz' simgenoHunksIII/SNPs.part.
ls simgenoHunksIII/SNPs.part* | cut -f2 -d"/" | cut -f3 -d"." > Phil.REAL.parts.SNP.txt 
ls simgenoHunksIII/HAPs.part* | cut -f2 -d"/" | cut -f3 -d"." > Phil.REAL.parts.HAP.txt

