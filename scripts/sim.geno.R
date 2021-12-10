
args = commandArgs(TRUE)
chr = as.character(args[1])   
Nsamples = as.numeric(args[2])
arrayID = as.numeric(args[3])
folder = as.character(args[4])

# chr = "Chr23"
# Nsamples=10
# arrayID=2
# folder="simgeno"

load(file.path(paste0("STI8_revII/", chr, "/RData"), paste0("EM.all.", chr, ".RData"))) ## or similar 
s = 1 ## assuming S=1 (normal way to run STITCH) 
sample_one_haplotype = function() { 
	n = rpois(n = 1, lambda = sum(-log(sigmaCurrent_m[, s]))) ## expected_number_of_recombinations 
	breaks = sort(sample(1:(nSNPs - 1), n, prob = -log10(sigmaCurrent_m) / sum(-log10(sigmaCurrent_m)))) 
	## between these, select hap dosages 
	starts = c(1, 1 + breaks) 
	ends = c(breaks, nSNPs) 
	prob = hapSumCurrent_tc[, starts, s] / N
	k = sapply(1:(length(breaks) + 1), function(i) sample(1:K, 1, prob=prob[,i]))
	SNPs = unlist(sapply(1:length(k),function(i) eHapsCurrent_tc[k[i], starts[i]:ends[i], s]))
	HAPs = unlist(sapply(1:length(k), function(i) rep(k[i],ends[i]-starts[i]+1)))
	list(SNPs,HAPs)
	} 
	
ff = function(x,Z){
	for(i in 1:(0.5*Nsamples)){
		Z[unlist(x[i]),i]<-1
		}
	Z
	}
	
blah=sample_one_haplotype()
goodSNP = info>0.40 & (estimatedAlleleFrequency > 0.01 & estimatedAlleleFrequency < 0.99)
NSNP=sum(goodSNP)
# generate samples
tempSNP = matrix(nrow=NSNP,ncol=Nsamples)
tempHAP = matrix(nrow=NSNP,ncol=Nsamples)
SS = seq(1,ncol(tempSNP),2)
for(i in 1:Nsamples){
	blah=sample_one_haplotype() 
	tempSNP[,i] = blah[[1]][goodSNP]
	tempHAP[,i] = blah[[2]][goodSNP]
	}

### SNPs
myID = (1:(0.5*Nsamples))+((arrayID-1)*0.5*Nsamples)
SNPgenos = round(tempSNP[,SS]+tempSNP[,SS+1],3)
myPOS=pos$POS[goodSNP]
if(arrayID==1){
	tempwrite=data.frame(CHR=rep(chr,nrow(SNPgenos)), POS=myPOS, SNPgenos)
	colnames(tempwrite)=c("CHR","POS",paste0("X",myID))
	}else{
	tempwrite=data.frame(SNPgenos)
	colnames(tempwrite)=paste0("X",myID)
	}
rm(tempSNP)
SNPfilename = paste0(folder,"/SNP_",arrayID,"_",chr,".txt")
write.table(tempwrite, SNPfilename, quote=FALSE, row.names = FALSE, col.names = TRUE)
rm (tempwrite)
gc()

### HAPs
# every 10th position
tempHAP = tempHAP[seq(1,NSNP,10),]
# row are haps
ZZ = matrix(rep(0,0.5*K*Nsamples),nrow=K)
bH1 = apply(tempHAP[,SS],1, ff, ZZ)
bH2 = apply(tempHAP[,SS+1],1, ff, ZZ)
rm (tempHAP)

# rows are individuals * positions (=SNPs), columns are 8 haplotypes predictors
temp = matrix(c(bH1),ncol=K,byrow=TRUE)+ matrix(c(bH2),ncol=K,byrow=TRUE)
rm (bH1)
rm (bH2)
temp = data.frame(temp)
colnames(temp) = paste("H",1:K,sep="")
temp2 = expand.grid(ID=myID,POS=pos$POS[goodSNP][seq(1,NSNP,10)])
temp2$CHR = rep(chr,nrow(temp2))
temp3 = data.frame(CHR=temp2$CHR, POS=temp2$POS, ID=temp2$ID, temp)
rm (temp)
rm (temp2)
HAPfilename = paste0(folder,"/HAP_",arrayID,"_",chr,".txt")

if(arrayID==1){
	write.table(temp3[order(temp3$ID, temp3$POS),],HAPfilename, quote=FALSE, row.names = FALSE, col.names = TRUE)
	}else{
	write.table(temp3[order(temp3$ID, temp3$POS),],HAPfilename, quote=FALSE, row.names = FALSE, col.names = FALSE)
	}

