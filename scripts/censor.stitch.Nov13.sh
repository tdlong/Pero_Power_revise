#!/bin/bash
#SBATCH --job-name=censor
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --array=1-12,14-24
#SBATCH --cpus-per-task=4  

KKK=(1 2 3 4 5 6 7 8a 8b 9 10 11 12 13 14 15 1621 17 18 19 20 22 23 X)
mychr=${KKK[$SLURM_ARRAY_TASK_ID-1]}
dir="STI8_revII/Chr"$mychr

module load samtools/1.10
module load bcftools/1.10.2

p=$dir"/stitch.Chr${mychr}.vcf.gz"
p2=$dir"/stitch.censor.Chr${mychr}.vcf.gz"
dir="PhilsampsII"

tabix $p
bcftools view -S Phillip.samples.txt -e 'INFO/INFO_SCORE<=0.4 || INFO/EAF < 0.01 || INFO/EAF > 0.99' -Oz -o $p2 $p 
tabix $p2
# tabix $p
bcftools view $p2 | bcftools query -f "%CHROM %POS[ %DS]\n" | gzip -c > $dir/SNPs.Chr${mychr}.censor.txt.gz

bcftools view $p2 | bcftools query -f "%CHROM %POS [%HD{0} %HD{1} %HD{2} %HD{3} %HD{4} %HD{5} %HD{6} %HD{7} ]\n" |\
	awk '(NR+5) % 10 == 0' | gzip -c > $dir/HAPs.Chr${mychr}.censor.txt.gz
zcat $dir/HAPs.Chr${mychr}.censor.txt | awk '{if(NR==FNR){a[$1]=$2;next}; for(i = 1; i <= 348; i++)\
	{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,a[i],$(8*(i-1)+3),$(8*(i-1)+4),$(8*(i-1)+5),\
	$(8*(i-1)+6),$(8*(i-1)+7),$(8*(i-1)+8),$(8*(i-1)+9),$(8*(i-1)+10))}}' samplenames.for.awk.txt - |\
	gzip -c >$dir/HAPs.Chr${mychr}.censor.rf.txt.gz

