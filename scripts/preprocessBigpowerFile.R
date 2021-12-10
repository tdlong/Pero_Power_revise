Spow = read.table("simgeno/bigpow_SNP.txt")
colnames(Spow) = c("CHR","PV","Gmod","CauseSNP","nHits","bestHit","whichBestHit","Nind")
library(tidyverse)
sSpow = Spow %>% 
	mutate(PV=recode(PV, PV10 = "PV1.0", PV25 = "PV2.5", PV50 = "PV5.0")) %>%
	group_by(CHR,PV,Gmod,Nind,CauseSNP) %>%	
	summarize(BH=max(bestHit), wBH = whichBestHit[which.max(bestHit)]) %>%
	group_by(CHR,PV,Gmod,Nind) %>%
	summarize(power=sum(BH>=7.5)/250, avgDistHit = mean((abs(wBH-CauseSNP))[BH>=7.5]))
sSpow$testFlav = rep("SNP",nrow(sSpow))

Hpow = read.table("simgeno/bigpow_HAP.txt")
colnames(Hpow) = c("CHR","PV","Gmod","CauseSNP","nHits","bestHit","whichBestHit","Nind")
library(tidyverse)
sHpow = Hpow %>% 
	mutate(PV=recode(PV, PV10 = "PV1.0", PV25 = "PV2.5", PV50 = "PV5.0")) %>%
	group_by(CHR,PV,Gmod,Nind,CauseSNP) %>%	
	summarize(BH=max(bestHit), wBH = whichBestHit[which.max(bestHit)]) %>%
	group_by(CHR,PV,Gmod,Nind) %>%
	summarize(power=sum(BH>=6.0)/250, avgDistHit = mean((abs(wBH-CauseSNP))[BH>=6.0]))
sHpow$testFlav = rep("HAP",nrow(sSpow))

write.table(rbind(sSpow,sHpow),"pow.pseudoInd.txt")


