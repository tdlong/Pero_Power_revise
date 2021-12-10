#!/bin/bash
#SBATCH --job-name=fq2bam
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --array=1-1275
#SBATCH --cpus-per-task=2   

module load bwa/0.7.8
module load samtools/1.10
module load bcftools/1.10.2

file="file.names.READ1.txt"
ref="/share/adl/tdlong/mouse_GWAS/data/ref/P.Leucopus.withUnreg.fa"
rawpath="/share/adl/tdlong/peromyscus/mouse_GWAS_raw_illumina/all_data"

# 4R086-L6-P01-TTCACATA-GTCATATT-READ1-Sequences.txt

read1=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1`
read2=`echo $read1 | sed 's/READ1/READ2/'`
readprefix=`echo $read1 | sed 's/-READ1-Sequences.txt//'`

bwa mem -t 2 -M $ref $rawpath/$read1 $rawpath/$read2 | samtools view -bS - > byplexBAM/$readprefix.bam
samtools sort byplexBAM/$readprefix.bam -o byplexBAM/$readprefix.sort.bam

