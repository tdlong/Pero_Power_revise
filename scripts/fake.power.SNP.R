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

sLOD = function(vv, yr){
	x = as.numeric(unlist(vv))
	N = length(x)
	YR = yr[1:N,]
	RSSfull = apply(resid(lm(YR ~ x))^2,2,sum)
	RSSnull = apply(resid(lm(YR ~ 1))^2,2,sum)
	F = ((RSSnull - RSSfull)) / (RSSfull/(N-2)) 
	LOD = -pf(F,1,N-2,lower.tail=FALSE,log.p=TRUE)/log(10)
	LOD
	}

yr = as.matrix(read.table("yr.2000.txt",header=TRUE))

snpin=paste0("simgenoHunks/SNPs.part.",part,".gz")
SS = fread(snpin,header=FALSE)    # CHR, POS, Ind1, Ind2 ... IndN
colnames(SS) = c("CHR","POS",paste0("X",1:2000))
SSs = SS %>% select(c(1:(Nind+2)))
rm(SS)
ptm = proc.time()
SS2 = SSs %>% group_by(CHR, POS) %>% nest() %>% mutate(LOD = map(data, sLOD, yr)) %>% select(-data)
rm(SSs)
ntests = nrow(SS2)
myLOD1 = SS2 %>% unnest(LOD)
rm(SS2)
myLOD1$test = rep(colnames(yr),ntests)
myLOD2 = myLOD1 %>% separate(test, c("CHRcause", "Gmod", "PV", "causeSNP"),sep="_") %>% 
	group_by(CHR, PV, Gmod, causeSNP) %>% summarize(nHIT = sum(LOD>7.5),bestHIT = max(LOD), whichBestHit = POS[which.max(LOD)])
myLOD2$Nind = rep(Nind,nrow(myLOD2))
cat("time to run SNP anova ",proc.time() - ptm,"\n")
tempfilename = paste0(myfolder,"/temppow_SNP_",Nind,"_",part,"_Chr23.txt")
write.table(myLOD2, tempfilename, quote=FALSE, row.names = FALSE, col.names = FALSE)

