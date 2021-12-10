#!/bin/bash
#SBATCH --job-name=fixbams
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --array=1-398
#SBATCH --cpus-per-task=2   

module load samtools/1.10
module load bamtools/2.5.1
module load bcftools/1.10.2
module load java/1.8.0
module load picard-tools/2.24.1
module load gatk/4.1.9.0

ref="/share/adl/tdlong/mouse_GWAS/data/ref/P.Leucopus.withUnreg.fa"
file="newIDs.txt"
ID=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1`

java -Xmx4G -jar /opt/apps/picard-tools/2.24.1/picard.jar AddOrReplaceReadGroups I=bysamBAM/$ID.bam O=fixBAM/$ID.RG.bam SORT_ORDER=coordinate RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=$ID RGSM=$ID VALIDATION_STRINGENCY=LENIENT
samtools index fixBAM/$ID.RG.bam
/opt/apps/gatk/4.1.9.0/gatk --java-options "-Xmx4G" MarkDuplicates -I fixBAM/$ID.RG.bam -O fixBAM/$ID.dedup.bam --REMOVE_DUPLICATES true -M bysamBAM/$ID.dedup.metrics.txt
samtools index fixBAM/$ID.dedup.bam
/opt/apps/gatk/4.1.9.0/gatk --java-options "-Xmx4G" HaplotypeCaller -R $ref -I fixBAM/$ID.dedup.bam --minimum-mapping-quality 30 -O fixBAM/$ID.dedup.vcf.gz -bamout fixBAM/$ID.IDR.bam
samtools index fixBAM/$ID.IDR.bam

