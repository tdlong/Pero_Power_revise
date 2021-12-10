#!/bin/bash
#SBATCH --job-name=stitch
#SBATCH -A tdlong_lab
#SBATCH -p standard          
#SBATCH --array=2-24
#SBATCH --cpus-per-task=32  

KKK=(1 2 3 4 5 6 7 8a 8b 9 10 11 12 13 14 15 1621 17 18 19 20 22 23 X)
mychr=${KKK[$SLURM_ARRAY_TASK_ID-1]}
dir="STI8_revII/Chr"$mychr
mkdir $dir

module load R/4.0.4
module load samtools/1.10

bamlist="bamlist2.Nov13.txt"
samnames="samplenames2.txt"
posfile="/share/adl/tdlong/mouse_GWAS/data/vcf/P.leu_"$mychr"_uniq.pos"
tempchr="Chr"$mychr
tdir="/share/adl/tdlong/peromyscus/mouse_GWAS_raw_illumina/makeBAM/TEMPDIR"
/share/adl/tdlong/peromyscus/mouse_GWAS_raw_illumina/makeBAM/STITCH/STITCH.R --chr=$tempchr --outputSNPBlockSize=5000 --tempdir=$tdir --bamlist=$bamlist --sampleNames_file=$samnames --posfile=$posfile --outputdir=$dir --K=8 --nGen=60 --nCores=16 --output_haplotype_dosages=TRUE

