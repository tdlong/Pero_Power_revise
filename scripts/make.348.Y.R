N = 348
yr = matrix(nrow=N,ncol=2*2*9*250)
fileOff = 0
lab=character(0)
BG = scan(file="BACKGROUND.348.txt",what=numeric(),sep=" ")
BG = BG[-length(BG)]
for(ff in c("Chr19","Chr23")){
	filename = paste0("PhilsampsII/SNPs.",ff,".censor.txt.gz")
	xx = read.table(filename,header=FALSE)     # CHR, POS, ... 2K columns of SNPs
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
		yr[,i + fileOff] =     sqrt(0.010)*scale(geno) + sqrt(0.990)*rnorm(N)
		yr[,i+(1*250) + fileOff] = sqrt(0.025)*scale(geno) + sqrt(0.975)*rnorm(N)
		yr[,i+(2*250) + fileOff] = sqrt(0.050)*scale(geno) + sqrt(0.950)*rnorm(N)
	# ten loci
		tempgeno_ten = apply(scale(t(xx[xx$POS %in% sample(tloc,min(10,length(tloc)),replace=FALSE),3:(N+2)])),1,sum)
		yr[,i+(3*250) + fileOff] = sqrt(0.010)*scale(tempgeno_ten) + sqrt(0.990)*rnorm(N)
		yr[,i+(4*250) + fileOff] = sqrt(0.025)*scale(tempgeno_ten) + sqrt(0.975)*rnorm(N)
		yr[,i+(5*250) + fileOff] = sqrt(0.050)*scale(tempgeno_ten) + sqrt(0.950)*rnorm(N)
	# all loci
		tempgeno_all = apply(t(xx[xx$POS %in% tloc,3:(N+2)]),1,sum)
		yr[,i+(6*250) + fileOff] = sqrt(0.010)*scale(tempgeno_all) + sqrt(0.990)*rnorm(N)
		yr[,i+(7*250) + fileOff] = sqrt(0.025)*scale(tempgeno_all) + sqrt(0.975)*rnorm(N)
		yr[,i+(8*250) + fileOff] = sqrt(0.050)*scale(tempgeno_all) + sqrt(0.950)*rnorm(N)
	# single locus + Background
		yr[,i+(9*250) + fileOff] =     sqrt(0.010)*scale(geno) + sqrt(0.049)*scale(BG) + sqrt(0.5)*rnorm(N)
		yr[,i+(10*250) +fileOff] = sqrt(0.025)*scale(geno) + sqrt(0.0475)*scale(BG) + sqrt(0.5)*rnorm(N)
		yr[,i+(11*250) +fileOff] = sqrt(0.050)*scale(geno) + sqrt(0.045)*scale(BG) + sqrt(0.5)*rnorm(N)
	# ten loci + Background
		tempgeno = apply(scale(t(xx[xx$POS %in% sample(tloc,min(10,length(tloc)),replace=FALSE),3:(N+2)])),1,sum)
		yr[,i+(12*250) + fileOff] = sqrt(0.010)*scale(tempgeno) + sqrt(0.049)*scale(BG) + sqrt(0.5)*rnorm(N)
		yr[,i+(13*250) + fileOff] = sqrt(0.025)*scale(tempgeno) + sqrt(0.0475)*scale(BG) + sqrt(0.5)*rnorm(N)
		yr[,i+(14*250) + fileOff] = sqrt(0.050)*scale(tempgeno) + sqrt(0.045)*scale(BG) + sqrt(0.5)*rnorm(N)
	# all loci + Background
		tempgeno = apply(t(xx[xx$POS %in% tloc,3:(N+2)]),1,sum)
		yr[,i+(15*250) + fileOff] = sqrt(0.010)*scale(tempgeno) + sqrt(0.049)*scale(BG) + sqrt(0.5)*rnorm(N)
		yr[,i+(16*250) + fileOff] = sqrt(0.025)*scale(tempgeno) + sqrt(0.0475)*scale(BG) + sqrt(0.5)*rnorm(N)
		yr[,i+(17*250) + fileOff] = sqrt(0.050)*scale(tempgeno) + sqrt(0.045)*scale(BG) + sqrt(0.5)*rnorm(N)
		}
	fileOff = fileOff + 2*9*250
	lab = c(lab,paste0(ff,"_SL_PV1.0_nB_",QTLloc),paste0(ff,"_SL_PV2.5_nB_",QTLloc),paste0(ff,"_SL_PV5.0_nB_",QTLloc),paste0(ff,"_TL_PV10_nB_",QTLloc),paste0(ff,"_TL_PV25_nB_",QTLloc),paste0(ff,"_TL_PV50_nB_",QTLloc),paste0(ff,"_AL_PV10_nB_",QTLloc),paste0(ff,"_AL_PV25_nB_",QTLloc),paste0(ff,"_AL_PV50_nB_",QTLloc))
	lab = c(lab,paste0(ff,"_SL_PV1.0_wB_",QTLloc),paste0(ff,"_SL_PV2.5_wB_",QTLloc),paste0(ff,"_SL_PV5.0_wB_",QTLloc),paste0(ff,"_TL_PV10_wB_",QTLloc),paste0(ff,"_TL_PV25_wB_",QTLloc),paste0(ff,"_TL_PV50_wB_",QTLloc),paste0(ff,"_AL_PV10_wB_",QTLloc),paste0(ff,"_AL_PV25_wB_",QTLloc),paste0(ff,"_AL_PV50_wB_",QTLloc))
	}	
yr=data.frame(yr)
sum(is.na(apply(yr,2,sum)))
colnames(yr) = lab
yr = round(yr,4)
fnout=paste0("yr.",N,".txt")
write.table(yr,file=fnout,quote=FALSE,row.names=FALSE)

