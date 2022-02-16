#Make rank abundance curves

#############
# the fungi #
#############

rm(list=ls())
library(tidyverse)
library(plyr)
library("wesanderson")

its_otus <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                     stringsAsFactors = F)

its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] -1
its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] / its_otus$ISD

samplingGroup <- gsub("X(\\d+_\\d+_\\d+)_.*", "\\1",its_otus$sample)

#Choose a random sample from each sampling group. 
chosenfew <- NA
k <- 1
for(i in samplingGroup){
  chosenfew[k] <- sample(which(samplingGroup == i),1)
  k <- k + 1
}

#split EN and ep
ens <- its_otus[grep("EN",its_otus$sample),]
eps <- its_otus[grep("EP",its_otus$sample),]

chosenfewTableEN <- ens[chosenfew,3:(length(ens)-2)] 
chosenfewTableEP <- eps[chosenfew,3:(length(eps)-2)] 

#repeat for bacteria

sixteens_otus <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                     stringsAsFactors = F)

sixteens_otus[,3:length(sixteens_otus)] <- sixteens_otus[,3:length(sixteens_otus)] -1
sixteens_otus[,3:length(sixteens_otus)] <- sixteens_otus[,3:length(sixteens_otus)] / sixteens_otus$ISD

samplingGroup <- gsub("X(\\d+_\\d+_\\d+)_.*", "\\1",sixteens_otus$sample)

#Choose a random sample from each sampling group. 
chosenfew <- NA
k <- 1
for(i in samplingGroup){
  chosenfew[k] <- sample(which(samplingGroup == i),1)
  k <- k + 1
}

#split EN and ep
ens16 <- sixteens_otus[grep("EN",sixteens_otus$sample),]
eps16 <- sixteens_otus[grep("EP",sixteens_otus$sample),]

chosenfewTableEN16 <- ens16[chosenfew,3:(length(ens16)-4)] 
chosenfewTableEP16 <- eps16[chosenfew,3:(length(eps16)-4)] 

pdf(width = 8, height = 8, file = "./visuals/rankabundanceCurvesISD.pdf")
par(mfrow = c(2,2))
hist(colSums(chosenfewTableEN, na.rm = T), 
     xlab = "abundance",
     main = "EN ITS")

hist(colSums(chosenfewTableEP, na.rm = T), 
     xlab = "abundance",
     main = "EP ITS")

hist(colSums(chosenfewTableEN16, na.rm = T), 
     xlab = "abundance",
     main = "EN 16S")

hist(colSums(chosenfewTableEP16, na.rm = T), 
     xlab = "abundance",
     main = "EP 16S")
dev.off()

#Do over for Hellinger transformed data. 


rm(list=ls())
library(tidyverse)
library(plyr)
library("wesanderson")

its_otus <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                     stringsAsFactors = F)

its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] -1
its_otus[,3:length(its_otus)]  <- vegan::decostand(its_otus[,3:(length(its_otus)-2)],
                                           method = "hellinger")

samplingGroup <- gsub("X(\\d+_\\d+_\\d+)_.*", "\\1",its_otus$sample)

#Choose a random sample from each sampling group. 
chosenfew <- NA
k <- 1
for(i in samplingGroup){
  chosenfew[k] <- sample(which(samplingGroup == i),1)
  k <- k + 1
}

#split EN and ep
ens <- its_otus[grep("EN",its_otus$sample),]
eps <- its_otus[grep("EP",its_otus$sample),]

chosenfewTableEN <- ens[chosenfew,3:(length(ens)-2)] 
chosenfewTableEP <- eps[chosenfew,3:(length(eps)-2)] 

#repeat for bacteria

sixteens_otus <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                          stringsAsFactors = F)

sixteens_otus[,3:length(sixteens_otus)] <- sixteens_otus[,3:length(sixteens_otus)] -1
sixteens_otus[,3:length(sixteens_otus)]  <- vegan::decostand(sixteens_otus[,3:(length(sixteens_otus)-4)],
                                                   method = "hellinger")
samplingGroup <- gsub("X(\\d+_\\d+_\\d+)_.*", "\\1",sixteens_otus$sample)

#Choose a random sample from each sampling group. 
chosenfew <- NA
k <- 1
for(i in samplingGroup){
  chosenfew[k] <- sample(which(samplingGroup == i),1)
  k <- k + 1
}

#split EN and ep
ens16 <- sixteens_otus[grep("EN",sixteens_otus$sample),]
eps16 <- sixteens_otus[grep("EP",sixteens_otus$sample),]

chosenfewTableEN16 <- ens16[chosenfew,3:(length(ens16)-4)] 
chosenfewTableEP16 <- eps16[chosenfew,3:(length(eps16)-4)] 

pdf(width = 8, height = 8, file = "./visuals/rankabundanceCurvesHella.pdf")
par(mfrow = c(2,2))
hist(colSums(chosenfewTableEN, na.rm = T), 
     xlab = "abundance",
     main = "EN ITS")

hist(colSums(chosenfewTableEP, na.rm = T), 
     xlab = "abundance",
     main = "EP ITS")

hist(colSums(chosenfewTableEN16, na.rm = T), 
     xlab = "abundance",
     main = "EN 16S")

hist(colSums(chosenfewTableEP16, na.rm = T), 
     xlab = "abundance",
     main = "EP 16S")
dev.off()
