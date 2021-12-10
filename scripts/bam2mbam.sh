#!/bin/bash
#SBATCH --job-name=bam2mbam
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --array=1-375
#SBATCH --cpus-per-task=2   

module load bwa/0.7.8
module load samtools/1.10
module load bamtools/2.5.1
module load bcftools/1.10.2
module load java/1.7
module load picard-tools/1.96

file="mouse.map.3.txt"
ID=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1 | cut -f1`
toPool=`head -n $SLURM_ARRAY_TASK_ID  $file | tail -n 1 | cut -f2`

bamtools merge $toPool -out bysamBAM/$ID.bam 
bamtools index -in bysamBAM/$ID.bam 

