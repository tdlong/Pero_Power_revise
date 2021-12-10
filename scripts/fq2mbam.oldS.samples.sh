#!/bin/bash
#SBATCH --job-name=fq2bam
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --array=1-48
#SBATCH --cpus-per-task=2   

module load bwa/0.7.8
module load samtools/1.10
module load bcftools/1.10.2
module load bamtools/2.5.1

file="file.names.oldS.cleanup.txt"
ref="/share/adl/tdlong/mouse_GWAS/data/ref/P.Leucopus.withUnreg.fa"
rawpath="/share/adl/tdlong/peromyscus/data"

IlluminaName=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1 | cut -f 1`
name=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1 | cut -f 2`

# bwa mem -t 2 -M $ref $rawpath/${IlluminaName}READ1-Sequences.txt.gz $rawpath/${IlluminaName}READ2-Sequences.txt.gz | samtools view -bS - > bysamBAM/$name.bam
# samtools sort bysamBAM/$name.bam -o byplexBAM/$name.sort.bam
# rm bysamBAM/$name.bam
# rm bysamBAM/$name.bam.bai
mv byplexBAM/$name.sort.bam bysamBAM/$name.bam
bamtools index -in bysamBAM/$name.bam 


