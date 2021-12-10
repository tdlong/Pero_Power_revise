library(tidyverse)
library(magrittr)
library(devtools)
library(lmtest)
# install_github("variani/lme4qtl")
library(lme4qtl)
# install_github("variani/matlm")
library(matlm)
# install_github("variani/wlm")
library(wlm)
library(Matrix)
library(data.table)

dofastscan = function(y, K, SNPtable){
	phen = data.frame(IND=row.names(SNPtable),Y=y)
	gfit_null_snp = lme4qtl::relmatLmer(formula = Y ~ (1 | IND), data = phen, relmat = list(IND = K), REML = TRUE)
	V = lme4qtl::varcov(gfit_null_snp, idvar = "IND")
	rm(gfit_null_snp)
	V[abs(V) < 1e-10] <- 0
	decomp = wlm::decompose_varcov(V, method = "evd", output = "all")
	rm(V)
	passoc_gls = matlm::matlm(formula = Y ~ 1, data = phen, pred = SNPtable, ids = IND, transform = decomp$transform, batch_size = 2000, verbose = 2)
	LODs = -log(passoc_gls$tab$pval) / log(10)                             # 0 values become NaN because log(0) = infinity
	rm(passoc_gls)
	LODs
	}

args = commandArgs(trailingOnly=TRUE)
myfolder = args[1]
part = args[2] 
# for testing
# part = "aa"
# myfolder="testtest"

K = read.table("kinship.popkin.txt")
cnames = as.character(row.names(K))	
row.names(K) = cnames
colnames(K) = cnames
K = as.matrix(K)

YR = read.table("yr.348.BT.txt",header=TRUE)
#   YR = YR[,1:25]
snpfile=paste0("simgenoHunksIII/SNPs.part.",part,".gz")
SS = fread(snpfile,header=FALSE)
colnames(SS) = c("CHR","POS",paste0("X",1:348)) 
SNPtable = t(SS[,-c(1,2)])
chrpos = SS[,c(1,2)]
rm(SS)
gc()
row.names(SNPtable) = cnames
temp = apply(YR, 2, dofastscan, K, SNPtable)
tempfilename1 = paste0(myfolder,"/temppow_LOD_",part,".txt")
tempfilename2 = paste0(myfolder,"/temppow_pow_",part,".txt")
write.table(cbind(chrpos,temp),tempfilename1,row.names=FALSE, quote=FALSE, col.names=FALSE)


