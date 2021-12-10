N = 2000
yr = matrix(nrow=N,ncol=9*250)
ff="Chr23"
filename = paste0("simgeno/SNP_Chr23.txt.gz")
xx = read.table(filename,header=TRUE)     # CHR, POS, ... 2K columns of SNPs
templabs = c("CHR","POS",paste0("X",1:N)) 
colnames(xx) = templabs
# check that each marker is polymorphic
tempsum=apply(xx[,3:(N+2)],1,sum)
xx=xx[tempsum>1 & tempsum<(2*N),]
filename2 = paste0("chooseSNP.",ff,".txt")
QTLloc = sort(scan(filename2,what=numeric()))

for(i in 1:250){
	qloc = QTLloc[i]
	tloc = xx$POS[xx$POS > (qloc - 50000) & xx$POS < (qloc + 50000)]
# single locus
	myelement = which.min(abs(xx$POS-qloc))    # this picks the closest SNP to qloc if qloc is monomorphic
	geno = c(unlist(xx[myelement,3:(N+2)]))
	yr[,i] =         sqrt(0.010)*scale(geno) + sqrt(0.990)*rnorm(N)
	yr[,i+(1*250)] = sqrt(0.025)*scale(geno) + sqrt(0.975)*rnorm(N)
	yr[,i+(2*250)] = sqrt(0.050)*scale(geno) + sqrt(0.950)*rnorm(N)
# ten loci
	tempgeno = apply(scale(t(xx[xx$POS %in% sample(tloc,min(10,length(tloc)),replace=FALSE),3:(N+2)])),1,sum)
	yr[,i+(3*250)] = sqrt(0.010)*scale(tempgeno) + sqrt(0.990)*rnorm(N)
	yr[,i+(4*250)] = sqrt(0.025)*scale(tempgeno) + sqrt(0.975)*rnorm(N)
	yr[,i+(5*250)] = sqrt(0.050)*scale(tempgeno) + sqrt(0.950)*rnorm(N)
# all loci
	tempgeno = apply(t(xx[xx$POS %in% tloc,3:(N+2)]),1,sum)
	yr[,i+(6*250)] = sqrt(0.010)*scale(tempgeno) + sqrt(0.990)*rnorm(N)
	yr[,i+(7*250)] = sqrt(0.025)*scale(tempgeno) + sqrt(0.975)*rnorm(N)
	yr[,i+(8*250)] = sqrt(0.050)*scale(tempgeno) + sqrt(0.950)*rnorm(N)
	}
lab = c(paste0(ff,"_SL_PV1.0_",QTLloc),paste0(ff,"_SL_PV2.5_",QTLloc),paste0(ff,"_SL_PV5.0_",QTLloc),paste0(ff,"_TL_PV10_",QTLloc),paste0(ff,"_TL_PV25_",QTLloc),paste0(ff,"_TL_PV50_",QTLloc),paste0(ff,"_AL_PV10_",QTLloc),paste0(ff,"_AL_PV25_",QTLloc),paste0(ff,"_AL_PV50_",QTLloc))
yr=data.frame(yr)
sum(is.na(apply(yr,2,sum)))
colnames(yr) = lab
yr = round(yr,4)
fnout=paste0("yr.",N,".txt")
write.table(yr,file=fnout,quote=FALSE,row.names=FALSE)

