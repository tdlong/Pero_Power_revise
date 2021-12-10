BT=read.table("k6.txt")
BT=BT[,c(8,9)]
cBT=colnames(BT)
BTbig = BT
for(i in 1:10){
	temp = BT[sample(1:348,348,replace=FALSE),]
	colnames(temp) = paste0(cBT,"_shuff_",i)
	BTbig = cbind(BTbig,temp)
	}
write.table(BTbig,"yr.348.BT.txt",quote=FALSE,row.names=FALSE)

K = read.table("kinship.popkin.txt")
PC15 = prcomp(K,retx=TRUE)$x[,1:15]
# read in phenotypes
YY = read.table("yr.348.BT.txt",header=TRUE)
temp = apply(YY,2,function(Y) lm(Y~PC15)$resid)
write.table(temp,"yr.348.BT.PC.txt")

