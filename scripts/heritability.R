library(lme4)
library(lme4qtl)

K = read.table("kinship.popkin.txt")
cnames = as.character(row.names(K))	
row.names(K) = cnames
colnames(K) = cnames
K = as.matrix(K)

YR = read.table("yr.348.BT.txt",header=TRUE)
phen = data.frame(Y = YR[,2], IND = cnames)
herBT = lme4qtl::relmatLmer(formula = Y ~ (1 | IND), data = phen, relmat = list(IND = K), REML = TRUE)
lme4qtl::VarProp(herBT)
# 6.5%

YR = read.table("yr.348.Weight.txt",header=TRUE)
phen = data.frame(Y = YR[,1], IND = cnames)
herW = lme4qtl::relmatLmer(formula = Y ~ (1 | IND), data = phen, relmat = list(IND = K), REML = TRUE)
lme4qtl::VarProp(herW)

