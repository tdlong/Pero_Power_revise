library(tidyverse)
library(data.table)
args = commandArgs(trailingOnly=TRUE)
myfolder = args[1]
Nind = as.numeric(args[2])      # number of individuals
part = args[3]                  # I will break the problem into 1Mb hunks

# for testing
# Nind = 100
# part = "aa"
# myfolder="simgeno"

mLOD = function(df, yr){
	drop_PC_less_than = 0.05
	x = df %>% select(starts_with("H")) %>% as.matrix()   # 8 columns
	pcs = stats::prcomp(x)
	pcs_filter <- pcs$x[, ((pcs$sdev^2 / sum(pcs$sdev^2)) > drop_PC_less_than)]
	YR = yr[1:nrow(pcs_filter),]	
	RSSfull = apply(resid(lm(YR ~ pcs_filter))^2,2,sum)
	RSSnull = apply(resid(lm(YR ~ 1))^2,2,sum)
	N = nrow(pcs_filter)
	p2 = ncol(pcs_filter)
	F = ((RSSnull - RSSfull)/p2) / (RSSfull/(N-p2-1)) 
	LOD = -pf(F,p2,N-p2-1,lower.tail=FALSE,log.p=TRUE)/log(10)
	LOD
	}

yr = as.matrix(read.table("yr.2000.txt",header=TRUE))
hapin=paste0("simgenoHunks/HAPs.part.",part,".gz")
HH = fread(hapin,header=FALSE)    # CHR, POS, ID, H1, ... H8
colnames(HH) = c("CHR","POS","ID",paste0("H",1:8))
HHs = HH %>% filter(ID<=Nind)
rm(HH)

ptm = proc.time()
HH2 = HHs %>% group_by(CHR, POS) %>% nest() %>% mutate(LOD = map(data, mLOD, yr)) %>% select(-data)
ntests = nrow(HH2)
myLOD1 = HH2 %>% unnest(LOD)
rm(HH2)
myLOD1$test = rep(colnames(yr),ntests)
myLOD2 = myLOD1 %>% separate(test, c("CHRcause", "Gmod", "PV", "causeSNP"),sep="_") %>%
	group_by(CHR, PV, Gmod, causeSNP) %>% summarize(nHIT = sum(LOD>6),bestHIT = max(LOD), whichBestHit = POS[which.max(LOD)])
myLOD2$Nind = rep(Nind,nrow(myLOD2))
cat("time to run HAP anova ",proc.time() - ptm,"\n")
tempfilename = paste0(myfolder,"/temppow_HAP_",Nind,"_",part,"_Chr23.txt")
write.table(myLOD2, tempfilename, quote=FALSE, row.names = FALSE, col.names = FALSE)
rm(myLOD1)
rm(myLOD2)
rm(HHs)
gc()

