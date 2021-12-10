library(popkin)
geno=read.table("hard.calls.txt",header=FALSE)
mynames=scan(file="Phillip.samples.txt",what=character())
K=popkin(as.matrix(geno))
colnames(K)=mynames
row.names(K)=mynames
pdf("kinship.pdf",height=5,width=5)
plot_popkin(K)
graphics.off()
# take top 15 PCs
PC15 = prcomp(K,retx=TRUE)$x[,1:15]
# read in phenotypes
YY = read.table("yr.348.txt",header=TRUE)
cc = colnames(YY)
YY2 = YY[,grepl("wB", cc)]
temp = apply(YY2,2,function(Y) lm(Y~PC15)$resid)
write.table(YY2,"yr.348.small.txt")
write.table(temp,"yr.348.PC.small.txt")
write.table(K,"kinship.popkin.txt")

