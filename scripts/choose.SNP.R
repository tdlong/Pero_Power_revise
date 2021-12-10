for(ff in c("Chr19","Chr23")){
	filename=paste0("PhilsampsII/SNPs.",ff,".censor.txt.gz")
	xx=read.table(filename,header=FALSE)
	templabs=c("CHR","POS",paste0("X",1:348))
	colnames(xx)=templabs
	# check each marker is polymorphic
	tempsum=apply(xx[,3:350],1,sum)
	xx=xx[tempsum>2 & tempsum < (2*348-2),]
	chooseSNP = sample(xx$POS,250,replace=FALSE)
	write(chooseSNP,file=paste0("chooseSNP.",ff,".txt"))
	}

