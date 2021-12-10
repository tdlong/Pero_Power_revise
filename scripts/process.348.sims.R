library(tidyverse)
library(data.table)

mmm = read.table("PhilscanSNP/tempbigpow_SNP.txt",header=FALSE)
colnames(mmm) = c("Pheno", "Chr_top_hit", "Pos_top_hit", "LOD_top_hit", "N_hit6")

myLOD = mmm %>%
	separate(Pheno, c("CauseChr", "Gmod", "PV", "BGV", "causeSNP"),sep="_") %>%
	mutate(causeSNP = as.numeric(causeSNP)) %>%
	mutate(PV=recode(PV, PV10 = "PV1.0", PV25 = "PV2.5", PV50 = "PV5.0")) %>%
	group_by(CauseChr, PV, Gmod, BGV, causeSNP) %>% 
	summarize(totalHits = sum(N_hit6), 
		bestHit = max(LOD_top_hit),
		whichBestChr = Chr_top_hit[which.max(LOD_top_hit)],
		whichBestHit = Pos_top_hit[which.max(LOD_top_hit)]) %>%
	mutate(correctChr = (whichBestChr==CauseChr),
		absDistKb = abs(whichBestHit - causeSNP)/1000)
# check that it is correct		
myLOD %>% ungroup() %>%
	group_by(CauseChr, PV, Gmod, BGV) %>% summarize(n())

write.table(myLOD,"PhilscanSNP/bigpow_SNP.txt")

myLOD2 = myLOD %>% ungroup() %>%
	group_by(CauseChr, PV, Gmod, BGV) %>%
	summarize(meanNHit = mean(totalHits[bestHit>7.5 & correctChr]),
		meanbestHit = mean(bestHit[bestHit>7.5 & correctChr]),
		percentHit=sum(bestHit>7.5)/250,
		percentHitCorrChr=sum(bestHit>7.5 & correctChr)/250,
		meandist=mean(absDistKb[bestHit>7.5 & correctChr]),
		mediandist=median(absDistKb[bestHit>7.5 & correctChr]))

myLOD2 %>% ungroup() %>% filter(CauseChr=="Chr19") %>% summarize(mean(percentHit))   # 3.6% so bang on

write.table(myLOD2,"PhilscanSNP/bigpow_SNP_summary.txt")
				
options("width"=200)
myLOD2 %>% print()

#with N=348 ...:
#	- power is low for PV = 2.5%
#	- power is low for PV = 5.0% but you would get a few
#	- localization is 1-3Mb, so not great, but OK for mapping with this N

sHdist %>% filter(causeCHR=="Chr19") %>% ungroup() %>% summarize(mean(percentHit)) = 1.78% so FPR is OK

correct CHR at:
	2.5% power = 9.2/3.6/4.0    for A/S/T genetic models
	5.0% power = 29.2/30.0/26.6 for A/S/T genetic models

I guess these are not significantly different between genetic models.  So I average ~6% @ 2.5% and 29% @ 5.0%

