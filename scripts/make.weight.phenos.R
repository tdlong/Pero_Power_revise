BT=read.table("k6.txt")
anova(lm(Weight~Age+Sex+timeDOB,data=BT))
anova(lm(log(Weight)~Sex+Age+Sex:Age,data=BT))
BT$RWeight = lm(log(Weight)~Sex+Age+Sex:Age,data=BT)$resid
quantile(BT$RWeight)
BTbig = matrix(nrow=nrow(BT),ncol=11)
BTbig[,1] = BT$RWeight
for(i in 2:11){
	BTbig[,i] = BTbig[sample(1:348,348,replace=FALSE),1]
	}
BTbig=data.frame(BTbig)
colnames(BTbig) = c("RWeight",paste0("RWeight_shuff_",1:10))
write.table(BTbig,"yr.348.Weight.txt",quote=FALSE,row.names=FALSE)
K = read.table("kinship.popkin.txt")
PC15 = prcomp(K,retx=TRUE)$x[,1:15]
# read in phenotypes
YY = read.table("yr.348.Weight.txt",header=TRUE)
temp = apply(YY,2,function(Y) lm(Y~PC15)$resid)
write.table(temp,"yr.348.Weight.PC.txt")

