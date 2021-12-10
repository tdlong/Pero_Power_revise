cd /share/adl/tdlong/peromyscus/mouse_GWAS_raw_illumina/makeBAM
# requires file.names.READ1.txt in helperfiles
sbatch fq2bam.sh
# requires mouse.map.3.txt in helperfiles
sbatch bam2mbam.sh
# required file.names.oldS.cleanup.txt in helperfiles
sbatch fq2mbam.oldS.samples.sh

# subsample S46 to ~1X
samtools view -bs 41.315 bysamBAM/S46.bam >temp.bam
samtools sort temp.bam -o bysamBAM/S46_1X.bam
mkdir fixBAM
# bam files bigger than 225000000 whic I estimate to be about 0.1X
ls bysamBAM/*.bam -alS | tr -s " " | awk '{if($5>225000000){print $9}}' | cut -d "/" -f2 | cut -d"." -f1 > newIDs.txt
# requires newIDs.txt in helperfiles
sbatch fixbam.sh

## run Stitch ....
ls fixBAM/*.IDR.bam | grep -v S46.IDR > bamlist2.txt
cat bamlist2.txt | sed 's:fixBAM/::' | sed 's:.IDR.bam::' > samplenames2.txt

tdir="/share/adl/tdlong/peromyscus/mouse_GWAS_raw_illumina/makeBAM/TEMPDIR"

mkdir STI8_revII
mkdir PhilsampsII
cat bamlist2.txt | sed 's/IDR/dedup/' >bamlist2.Nov13.txt
sbatch stitch.Nov13.sh

### identify samples with Bleeding time data
p="STI8_revII/Chr23/stitch.Chr23.vcf.gz"
zcat $p | head -n 12 | tail -n 1 | cut -f 10- >geno.samples.txt

# required k5.txt in helperfiles
R
k5=read.table("k5.txt",header=TRUE)
g=scan(file="geno.samples.txt",what=character(),sep="\t")
gtemp = g[g %in% as.character(k5$ID)]
k6 = k5[match(gtemp,as.character(k5$ID)),]
write.table(k6,"k6.txt")
for(i in 1:length(gtemp)){
	cat(gtemp2[i],"\n",file=Phillip.samples.txt,append=TRUE)
	}

wc -l Phillip.samples.txt # 348
cat Phillip.samples.txt | tr -d ' ' > temp.txt
mv temp.txt Phillip.samples.txt

# this will create new stitch files that are censored on MAF >1% and info>0.4
# it also creates simpler snp and hap tables (gzipped) in Philsamps folder
p2="STI8_revII/Chr17/stitch.Chr17.vcf.gz"
bcftools view -h $p2 | tail -n 1 | cut -f 10- | tr "\t" "\n" | awk '{printf("%s\t%s\n",NR,$1)}' >samplenames.for.awk.txt
sbatch censor.stitch.Nov13.sh


# Check RNAseq versus imputation calls
# RNAseq = spleen.genos.txt in helperfiles
# Chr	Pos	R20949	R20819	R20733	R20495	R20483	R20695	R20688	R20807	R20702	R20694	R20732
# Chr10	33534	0	0	1	0	1	0	1	0	0	1	0

cat helperfiles/spleen.genos.txt | grep Chr23 | cut -f2,3 | awk '{printf("%s\t%s\n",$1,$2)}' | sort -k1,1 > RNA_Chr23_20949.txt
zcat PhilsampsII/SNPs.Chr23.censor.txt.gz | cut -d " " -f2,83,84,85,86,87 | awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$5,$6)}' | sort -k1,1 > impute_Chr23_20949.txt
join -j 1 RNA_Chr23_20949.txt impute_Chr23_20949.txt > join_20949.txt 


########### PSEUDO INDIVIDUALS ###########

# simulate genotype and consolidate
# simulated samples are numbered 1-2000
mkdir simgeno
sbatch sim.geno.sh
sbatch consolidate.sh
ls -alth simgeno
# Nov 14 10:56 SNP_Chr23.txt.gz
# Nov 14 10:53 HAP_Chr23.txt.gz

# choose SNPs to be QTL
sbatch choose.SNP.sh
ls -alth chooseSNP*
# Nov 14 11:02 chooseSNP.Chr23.txt
# Nov 14 11:01 chooseSNP.Chr19.txt

### MAKE PHENOTYPES
### only Chr23 only no Background.
Rscript make.fake.Y.R

# check
cat yr.2000.txt | head -n5 | cut -f1-5 -d" "
Chr23_SL_PV1.0_619253 Chr23_SL_PV1.0_741131 Chr23_SL_PV1.0_1461090 Chr23_SL_PV1.0_1504176 Chr23_SL_PV1.0_1740882
-0.1331 -0.7365 -1.8932 1.1247 0.8527
0.5145 0.3026 -1.1191 -1.5115 0.7723
0.4417 1.1451 -0.8034 -0.4789 -1.4628
0.7772 0.9878 -1.5336 0.4282 0.5167

## break Chr23 into 8 pieces, write a file to break the scan into parts
mkdir simgenoHunks
zcat simgeno/SNP_Chr23.txt.gz | wc -l    # 352578
zcat simgeno/HAP_Chr23.txt.gz | wc -l    # 70516001
sbatch break.chr23.hunks.sh
ls simgenoHunks/SNPs.part* | cut -f2 -d"/" | cut -f3 -d"." > temp.parts.SNP.txt
ls simgenoHunks/HAPs.part* | cut -f2 -d"/" | cut -f3 -d"." > temp.parts.HAP.txt
# eg
# simgenoHunks/HAPs.part.aa.gz 
# simgenoHunks/SNPs.part.aa.gz 

R
Nind=c(100,250,348,500,750,1000,1500,2000)
part=scan("temp.parts.SNP.txt",what=character())
temp=expand.grid(Nind,part)
write.table(temp,"Chr23.hunks.SNP.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
part=scan("temp.parts.HAP.txt",what=character())
temp=expand.grid(Nind,part)
write.table(temp,"Chr23.hunks.HAP.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

wc -l Chr23.hunks.HAP.txt #320
wc -l Chr23.hunks.SNP.txt #288

## some jobs may need to be rerun with more memory, if you have infinite resources just give the jobs more memory
sbatch fake.geno.HAP.sh   # 1 core is enough, these are fast maybe < 30 min
sbatch fake.geno.SNP.sh   # 2 core is enough, these can take some time for bigger N, 9-10 hours per job

cat simgeno/temppow_HAP* > simgeno/bigpow_HAP.txt
rm simgeno/temppow_HAP*
cat simgeno/temppow_SNP* > simgeno/bigpow_SNP.txt
rm simgeno/temppow_SNP*

sbatch preprocessBigpowerFile.sh
# output = pow.pseudoInd.txt 

Rscript plot.pseudo.power.R
# output = pseudo.ind.pdf
# and some other information, look at the script

################################################
################################################
##############   348 samples ###################
################################################
################################################

########### 348 INDIVIDUALS ###########

# make background
# sum of every 400th SNP (starting at 200)
sbatch make.background.sh
dir="PhilsampsII"
cat $dir/SNPs.Chr*.BG.txt | awk '{for(i=1; i<=NF; i++) total[i] += $i}; END {for(i=1; i<=NF; i++){printf "%s ", total[i]};print ""}' > BACKGROUND.348.txt
rm $dir/SNPs.*.BG.txt 

### MAKE PHENOTYPES
### the phenotypes for the real samples have 250 replicates, BG vs. noBG, 3 levels of %var, Chr19 vs. Chr23 causative
### X does not contribute
Rscript make.348.Y.R

# check
cat yr.348.txt | head -n5 | cut -f1-5 -d" "
Chr19_SL_PV1.0_nB_1692967 Chr19_SL_PV1.0_nB_1726527 Chr19_SL_PV1.0_nB_2333890 Chr19_SL_PV1.0_nB_2419764 Chr19_SL_PV1.0_nB_2426818
0.9231 -0.4701 0.4629 0.7002 0.6465
-1.1028 -1.5988 1.3944 0.51 0.4264
-2.2948 -0.6877 -0.8371 -1.1519 -0.8354
-0.3833 -1.6754 0.3919 -0.7429 0.6542

#  1. format of yr.348.txt & Phillip.samples.txt
#	- columns are different sets of 348 Y's
#	- rows are different individuals, unlabeled, but given by Phillip.samples.txt
#	- col headers describe the condition, eg. "Chr19_SL_PV1.0_nB_1008824" are phenotypes for
#		-cause SNP at Chr19:1008825
#		-the SL genetic model (S,A,T are single, all, ten causative SNPs)
#		-PV1.0 means QTL accounts for 1% of total variation (also 2.5 & 5.0 tested)
#		-nB tells me n0 background (only SNP and Ve), compared to wB where 50% of Vt is genetic
###  Real phenotypes have to be in the same order as SNPs

### break the genotypes into equal sized hunks
### also extract SNPs to make kinship matrix

mkdir simgenoHunksII
sbatch breakgenome.sh

### output
PhilsampsII/HAP.all.sort.txt  # concatenate genome
PhilsampsII/SNP.all.sort.txt
simgenoHunksII/*              # above in pieces
Phil.parts.SNP.txt   #list of 322 indicies for subsetting
Phil.parts.HAP.txt   #322
hard.calls.txt       # every 631st SNP as hardcall to make kinship matrix

# number of SNPs post filtering
wc -l PhilsampsII/SNP.all.sort.txt   # -> 16087368

sbatch kinPCsmall.sh

### output
kinship.pdf          # kinship plot
yr.348.small.txt     # phenotypes for 348 scan (all wB)
yr.348.PC.small.txt  # 1st 15 PCs using kinship
kinship.popkin.txt   # kinship matrix

mkdir PhilscanSNP
mkdir PhilscanHAP

sbatch HAPer.II.sh
cat PhilscanHAP/temppow_HAP* | gzip -c > PhilscanHAP/bigpow_HAP.txt.gz
rm PhilscanHAP/temppow_HAP*
cat PhilscanHAP/tempLOD_HAP* | gzip -c > PhilscanHAP/bigLOD_HAP.txt.gz   

R
library(data.table)
library(tidyverse)
Hpow = fread("PhilscanHAP/bigpow_HAP.txt.gz",header=FALSE)
colnames(Hpow) = c("testCHR","causeCHR","PV","Gmod","BG","CauseSNP","nHits","bestHit","whichBestHit")
sHdist = Hpow %>% mutate(PV=recode(PV, PV10 = "PV1.0", PV25 = "PV2.5", PV50 = "PV5.0")) %>% 
	group_by(causeCHR,CauseSNP,PV,Gmod,BG) %>%
	slice(which.max(bestHit)) %>% filter(nHits>=1) %>%
	mutate(correctCHR = (testCHR==causeCHR)) %>%
	mutate(distFromCause = abs(whichBestHit-CauseSNP)) %>%
	group_by(causeCHR,PV,Gmod,BG) %>%
	summarize(percentHit=n()/250,percentCorrChr=sum(correctCHR)/250,mdist=mean(distFromCause[correctCHR]))
sHdist$testFlav = rep("HAP",nrow(sHdist))
write.table(sHdist,"temp.pow.348.HAPonly.txt")


## SNPs were harder as the file sizes are huge
## this has got it
## but now the final step is separated from generating the LOD files
## the first script generates the LODs using 2 cores and 12 hours
## the second (process) script uses xx cores but only xx minutes

sbatch SNPer.sim4500.sh
## gzip the LODs ... I submitted this as job as it could take days...
cat PhilscanSNP/temppow_LOD_* | gzip -c > PhilscanSNP/bigLOD_SNP.txt.gz

## partially process raw LODs to extract some information
## occasionally python is preferred to R!
head -n 1 yr.348.small.txt > colnames.yr.348.small.txt
sbatch SNPer.sim4500.process.sh

# in = PhilscanSNP/temppow_LOD_aa.txt   -> Chr, Pos, Pheno1, Pheno2,..., Pheno4500
# out =  PhilscanSNP/temppow_pow_aa.txt -> Pheno, Chr_top_hit, Pos_top_hit, LOD_top_hit, N_hit_g7.5  with 4500 rows = phenotypes
cat PhilscanSNP/temppow_pow* > PhilscanSNP/tempbigpow_SNP.txt
Rscript process.348.sims.R

# writing to PhilscanSNP and PhilscanHAP
# bigpow are summaries
# bigLOD are per marker LOD scores
# note the formating of the LOD files
# bigLOD_SNP.txt.gz            # all the LOD scores for SNPs {CHR,POS,LOD1,LOD2,...} ... head -n1 yr.348.small.txt gives names of columns 3..ncol
# bigLOD_HAP.txt.gz            # all the LOD scores for HAPs {CHR,POS,LOD,PHENO_FROM_WHERE}

# summaries =
mv temp.pow.348.HAPonly.txt PhilscanHAP/bigpow_HAP_summary.txt
PhilscanSNP/bigpow_SNP.txt
PhilscanSNP/bigpow_SNP_summary.txt

###############################
### REAL PHENOTYPES
###############################
### bleeding time
###############################
Rscript make.BT.phenos.R
# output = yr.348.BT.txt & yr.348.BT.PC.txt

# add Chr 19 back to the scan
mkdir simgenoHunksIII
sbatch breakgenome.REAL.sh

### output
simgenoHunksIII/HAP.all.sort.txt  # concatenate genome
simgenoHunksIII/SNP.all.sort.txt
simgenoHunksIII/*              # above in pieces
Phil.REAL.parts.HAP.txt   #list of 338 indicies for subsetting
Phil.REAL.parts.SNP.txt   #338

mkdir PhilRealSNP
mkdir PhilRealHAP

sbatch HAPer.BT.sh 
sbatch SNPer.BT.sh

# writing to PhilscanSNP and PhilscanHAP
echo "CHR POS LOD dataset" > PhilRealHAP/bigLOD_HAP.txt
cat PhilRealHAP/tempLOD_HAP_* >>  PhilRealHAP/bigLOD_HAP.txt
rm PhilRealHAP/tempLOD_HAP_*

head -n1 yr.348.BT.txt
echo "CHR POS normBleed residNormBleed normBleed_shuff_1 residNormBleed_shuff_1 normBleed_shuff_2 residNormBleed_shuff_2 normBleed_shuff_3 residNormBleed_shuff_3 normBleed_shuff_4 residNormBleed_shuff_4 normBleed_shuff_5 residNormBleed_shuff_5 normBleed_shuff_6 residNormBleed_shuff_6 normBleed_shuff_7 residNormBleed_shuff_7 normBleed_shuff_8 residNormBleed_shuff_8 normBleed_shuff_9 residNormBleed_shuff_9 normBleed_shuff_10 residNormBleed_shuff_10" > PhilRealSNP/bigLOD_SNP.txt
cat PhilRealSNP/temppow_LOD_* >> PhilRealSNP/bigLOD_SNP.txt
rm PhilRealSNP/temppow_LOD_*

R
library(tidyverse)
library(data.table)
temp=fread("PhilRealSNP/bigpow_SNP.txt",header=TRUE)
out = apply(temp[,3:ncol(temp)],2,function(x) sum(x > 7.5))
data.frame(row.names=names(out),hits=out)
temp2=fread("PhilRealHAP/bigpow_HAP.txt",header=TRUE)
out2 = temp2 %>% group_by(dataset) %>% summarize(nhit = sum(LOD > 6.0))

###############################
### log(weight)
###############################

Rscript make.weight.phenos.R

mkdir PhilWeightSNP
mkdir PhilWeightHAP
sbatch HAPer.Weight.sh
sbatch SNPer.Weight.sh

# writing to PhilscanSNP and PhilscanHAP
echo "CHR POS LOD dataset" > PhilWeightHAP/bigLOD_HAP.txt
cat PhilWeightHAP/tempLOD_HAP_* >> PhilWeightHAP/bigLOD_HAP.txt
rm PhilWeightHAP/tempLOD_HAP_*
head -n1 yr.348.Weight.txt
echo "CHR POS RWeight RWeight_shuff_1 RWeight_shuff_2 RWeight_shuff_3 RWeight_shuff_4 RWeight_shuff_5 RWeight_shuff_6 RWeight_shuff_7 RWeight_shuff_8 RWeight_shuff_9 RWeight_shuff_10" > PhilWeightSNP/bigLOD_SNP.txt
cat PhilWeightSNP/temppow_LOD_* >> PhilWeightSNP/bigLOD_SNP.txt
rm PhilWeightSNP/temppow_LOD_*

R
library(tidyverse)
library(data.table)
temp=fread("PhilWeightSNP/bigLOD_SNP.txt",header=TRUE)
out = apply(temp[,3:ncol(temp)],2,function(x) sum(x > 7.5))
out = apply(temp[,3:ncol(temp)],2,function(x) sum(x > 6.0))
data.frame(row.names=names(out),hits=out)

temp2=fread("PhilWeightHAP/bigLOD_HAP.txt",header=TRUE)
out2 = temp2 %>% group_by(dataset) %>% summarize(nhit = sum(LOD > 6.0))
out2

temp2 %>% filter(LOD>6)

CHR      POS      LOD dataset
Chr3 69841667 6.214956 RWeight = NC_051065.1:69841667

#######  heritability

Rscript heritability.R

####################
# Where is stuff
####################
# SNPcalls  (these are new filtered calls)
# Chr pos call#1, #2, etc
$dir/PhilsampsII/SNPs.Chr*.censor.txt.gz
$dir/PhilsampsII/HAPs.Chr*.censor.rf.txt.gz 
# the sample names
$dir/Phillip.samples.txt
# the phenotypes
$dir/k6.txt
# simulated phenotypes for the 348 individuals (interesting cases could be pulled out and run for Chr23)
#  format of yr.348.txt & Phillip.samples.txt
#	- columns are different sets of 348 Y's
#	- rows are different individuals, unlabels, but given by Phillip.samples.txt
#	- col headers describe the condition, eg. "Chr19_SL_PV1.0_nB_1008824" are phenotypes for
#		-cause SNP at Chr19:1008825
#		-the SL genetic model (S,A,T are single, all, ten causative SNPs)
#		-PV1.0 means QTL accounts for 1% of total variation (also 2.5 & 5.0 tested)
#		-nB tells me n0 background (only SNP and Ve), compared to wB where 50% of Vt is genetic
$dir/yr.384.small.txt               # 348 individuals, all wB, 4500 combinations (2 chr * 3 gen model * 3 %var * 250 replicate loci)
$dir/yr.348.PC.small.txt            # same as above but 1st 15 PCs based on kinship removed
$dir/yr.2000.txt                    # simulated phenotypes for pseudo individuals scan
$dir/yr.348.BT.txt                  # 348 individuals, norm bleeding time or resid norm bleedtime + 10 perms of same
$dir/yr.348.BT.PC.txt               # same as above but 1st 15 PCs based on kinship removed
$dir/yr.348.Weight.txt              # 348 weight
$dir/yr.348.Weight.PC.txt           # 348 weight PC correcrted

# results of my scan of the 4500 columns above that are "wB"
# chr, pos, pheno#1, pheno#2, etc.
# summaries
# columns are: c("testCHR","causeCHR","PV","Gmod","BG","CauseSNP","nHits","bestHit","whichBestHit")
$dir/PhilscanHAP/bigpow_HAP.txt.gz
$dir/PhilscanSNP/bigpow_SNP.txt.gz
# all LODs ...
# they do not have headers 
$dir/PhilscanSNP/bigLOD_SNP.txt.gz            # all the LOD scores for SNPs {CHR,POS,LOD1,LOD2,...} ... head -n1 yr.348.small.txt gives names of columns 3..ncol
$dir/PhilscanHAP/bigLOD_HAP.txt.gz            # all the LOD scores for HAPs {CHR,POS,LOD,PHENO_FROM_WHERE}

# the kinship matrix made with SNPs
$dir/kinship.popkin.txt
# results of scans of BT or permuted BT, all positions all LODs 
# residNormBT is the main character, it can be compared to permutations for HAP or SNP
# qqplot or Manhattan plots
PhilRealSNP/bigLOD_SNP.txt      # CHR, POS, LOD_A, LOD_B, etc   ... head -n1 yr.348.BT.txt gives names of columns 3..ncol
PhilRealHAP/bigLOD_HAP.txt      # CHR, POS, LOD, dataset  ... so more tidy

# platelet GO category P.leucopus hits
helperfiles/PlateletHits.txt
# why does SCGB not use natural names...converter table
helperfiles/SCBG2normal.chromosomeName.txt
# a less terse (more impenetrable) set of notes
helperfiles/muchMoreExtensiveNotes.txt
