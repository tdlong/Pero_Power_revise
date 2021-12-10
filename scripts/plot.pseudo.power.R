library(tidyverse)
pow = read.table("pow.pseudoInd.txt")
library(grid)
library(gridExtra)

A = ggplot(pow %>% filter(Gmod=="SL"), aes(x = Nind, y = power)) + 
  geom_line(aes(color = PV, linetype = testFlav)) + 
  labs(y="power", x = "N individual") +
  theme_gray(base_size = 9)
B = ggplot(pow %>% filter(Gmod=="TL"), aes(x = Nind, y = power)) + 
  geom_line(aes(color = PV, linetype = testFlav)) + 
  labs(y="power", x = "N individual") +
  theme_gray(base_size = 9)
C = ggplot(pow %>% filter(Gmod=="AL"), aes(x = Nind, y = power)) + 
  geom_line(aes(color = PV, linetype = testFlav)) + 
  labs(y="power", x = "N individual") +
  theme_gray(base_size = 9)

# dotted = SNP, solid = HAP
# top to botton SL, TL, AL
pdf("pseudo.ind.pdf",height=5,width=4)
grid.arrange(A+theme(legend.position='hidden'), B+theme(legend.position='hidden'),
             C+theme(legend.position='hidden'),
             layout_matrix=matrix(c(1,2,3)),ncol=1)
graphics.off()
       
pow %>% filter(Nind==348) %>% select(PV,Gmod,testFlav,power)

PV Gmod testFlav power
PV1.0   AL      SNP 0.000
PV1.0   SL      SNP 0.000
PV1.0   TL      SNP 0.000
PV2.5   AL      SNP 0.000
PV2.5   SL      SNP 0.000
PV2.5   TL      SNP 0.008
PV5.0   AL      SNP 0.156
PV5.0   SL      SNP 0.164
PV5.0   TL      SNP 0.112
PV1.0   AL      HAP 0.000
PV1.0   SL      HAP 0.004
PV1.0   TL      HAP 0.008
PV2.5   AL      HAP 0.008
PV2.5   SL      HAP 0.012
PV2.5   TL      HAP 0.008
PV5.0   AL      HAP 0.108
PV5.0   SL      HAP 0.140
PV5.0   TL      HAP 0.140

pow$avD = pow$avgDistHit/1000
pow %>% filter(PV == "PV5.0"  & Gmod != "SL" & Nind >= 348) %>%
	select(-c(CHR,PV,power,avgDistHit)) %>%
	pivot_wider(names_from = testFlav, values_from = avD)

Gmod   Nind   SNP   HAP
 AL      348 626.   72.9
 AL      500 100.   60.2
 AL      750  45.6  52.9
 AL     1000  44.6  47.4
 AL     1500  37.3  37.0
 AL     2000  34.1  32.3
 TL      348 781.  169. 
 TL      500 379.   63.9
 TL      750  98.9  64.7
 TL     1000  49.8  56.4
 TL     1500  76.7  46.3
 TL     2000  72.5  73.4
 
#Gmod = AL (all) or TL (ten cause).  I ignore SL as localization can be extremely good.
#SNP = average distance in kb between the most significant marker and the causative locus (conditional on a hit) for SNP-based tests
#HAP = average distance in kb between the most significant marker and the causative locus (conditional on a hit) for HAP-based tests

#...when you get a hit, localization can be very good for 5% QTL, even when N=348

pow %>% filter(PV == "PV2.5"  & Gmod != "SL" & Nind >= 348) %>%
	select(-c(CHR,PV,power,avgDistHit)) %>%
	pivot_wider(names_from = testFlav, values_from = avD)

# The same trend hold for 2.5% var ... provided samples sizes are larger (say N >=1000)

Gmod   Nind    SNP     HAP
AL      348 3998.  21012. 
AL      500 2522.   3812. 
AL      750  366.   1249. 
AL     1000  270.     67.8
AL     1500   47.6   138. 
AL     2000   45.9    50.9
TL      348  909.     29.2
TL      500 1034.   2075. 
TL      750  441.     60.8
TL     1000  312.     77.1
TL     1500   57.0    62.0
TL     2000   91.3    56.7
