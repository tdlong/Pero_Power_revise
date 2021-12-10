library(tidyverse)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
myfolder = args[1]
part = args[2] 
# for testing
# part = "aa"
# myfolder="testtest"


mLOD = function(df, yr){
	drop_PC_less_than = 0.05
	x = df %>% select(starts_with("H")) %>% as.matrix()   # 8 columns
	pcs = stats::prcomp(x)
	pcs_filter <- pcs$x[, ((pcs$sdev^2 / sum(pcs$sdev^2)) > drop_PC_less_than)]
	RSSfull = apply(resid(lm(yr ~ pcs_filter))^2,2,sum)
	RSSnull = apply(resid(lm(yr ~ 1))^2,2,sum)
	N = nrow(as.matrix(pcs_filter))
	p2 = ncol(as.matrix(pcs_filter))
	F = ((RSSnull - RSSfull)/p2) / (RSSfull/(N-p2-1)) 
	LOD = -pf(F,p2,N-p2-1,lower.tail=FALSE,log.p=TRUE)/log(10)
	LOD
	}

yr = as.matrix(read.table("yr.348.PC.small.txt",header=TRUE))
hapfile=snpfile=paste0("simgenoHunksII/HAPs.part.",part,".gz")
HH = fread(hapfile,header=FALSE)
colnames(HH) = c("CHR","POS","ID","H1","H2","H3","H4","H5","H6","H7","H8") 

HH2 = HH %>% group_by(CHR, POS) %>% nest() %>% mutate(LOD = map(data, mLOD, yr)) %>% select(-data)
ntests = nrow(HH2)
myLOD1 = HH2 %>% unnest(LOD)
rm(HH2)
rm(HH)
myLOD1$test = rep(colnames(yr),ntests)
myLOD2 = myLOD1 %>% separate(test, c("CauseChr", "Gmod", "PV", "BGV", "causeSNP"),sep="_") %>%
	group_by(CHR, CauseChr, PV, Gmod, BGV, causeSNP) %>% summarize(nHIT = sum(LOD>6),bestHIT = max(LOD), whichBestHit = POS[which.max(LOD)])
tempfilename1 = paste0(myfolder,"/tempLOD_HAP_",part,".txt")
# CHR, POS, pheno_replicate = [c("CauseChr", "Gmod", "PV", "BGV", "causeSNP"),sep="_"], LOD
tempfilename2 = paste0(myfolder,"/temppow_HAP_",part,".txt")
write.table(myLOD1, tempfilename1, quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(myLOD2, tempfilename2, quote=FALSE, row.names = FALSE, col.names = FALSE)

