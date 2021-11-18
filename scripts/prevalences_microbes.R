#############
# the fungi #
#############

rm(list=ls())
library(tidyverse)
library(plyr)
library("wesanderson")

#SINTAX
its <- read.table("processedData/smallmem97_ITS.sintax", fill = T)
head(its)

its_otus <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                     stringsAsFactors = F)
length(names(its_otus)) -4

its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] -1
its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] / its_otus$ISD

its_meta <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                     stringsAsFactors = F)

its_otus$sample <- gsub("X", "", its_otus$sample)

merged_dat <- merge(its_meta, its_otus, by.x = "sample", by.y = "sample", all.y = T)
dim(merged_dat)
dim(its)
dim(its_otus)
names(merged_dat)[200:205]

#recode the names of merged dat for each Zotu as the taxonomy hypothesis
for(i in 204:length(merged_dat)){
  if(any(grepl(names(merged_dat)[i], as.character(its$V1)))){
    names(merged_dat)[i] <- paste(names(merged_dat)[i],
                                  its$V4[its$V1 == names(merged_dat)[i]],
                                  sep = "_")
  }
}

#binarize
merged_dat[,grep("Zotu",names(merged_dat))] <- ifelse(merged_dat[,grep("Zotu",names(merged_dat))] > 0, 1, 0)

#agg by host 
hits <- list()
k <- 1
for(i in grep("Zotu",names(merged_dat))){
  hits[[k]] <- aggregate(merged_dat[,i] ~ merged_dat$taxon_final, FUN = sum)
  k <- k +1
}

#Overall prevalence
prevalenceF <- colSums(merged_dat[,grep("Zotu",names(merged_dat))], na.rm = T)
summary(prevalenceF)
hist(prevalence)

#combine rows by host, rebinarize, and redetermine prevalence to see which things are in the most hosts
counts <- merged_dat[,grep("Zotu",names(merged_dat))]

agg_hostEN <- data.frame(matrix(nrow = length(unique(merged_dat$taxon_final)),
                              ncol = 1))
k <- 1
for(i in 1:length(counts)){
  agg_hostEN[,k] <- aggregate(counts[merged_dat$compartment == "EN",i] ~ merged_dat$taxon_final[merged_dat$compartment == "EN"], 
                            FUN = sum)[,2]
  colnames(agg_hostEN)[k] <- names(counts)[i]
  k <- k + 1
}
agg_hostEN <- data.frame(unique(merged_dat$taxon_final),
                       agg_hostEN)

agg_hostEN[,2:length(agg_hostEN)] <- ifelse(agg_hostEN[,2:length(agg_hostEN)] > 0, 1, 0)
prevHostEN <- colSums(agg_hostEN[,2:length(agg_hostEN)])
summary(prevHostEN)
hist(prevHostEN)
dim(agg_hostEN)
table(prevHostEN > 50)

prevHostEN[prevHostEN > 50]

#FOr EP
agg_hostEP <- data.frame(matrix(nrow = length(unique(merged_dat$taxon_final)),
                                ncol = 1))
k <- 1
for(i in 1:length(counts)){
  agg_hostEP[,k] <- aggregate(counts[merged_dat$compartment == "EP",i] ~ merged_dat$taxon_final[merged_dat$compartment == "EP"], FUN = sum)[,2]
  colnames(agg_hostEP)[k] <- names(counts)[i]
  k <- k + 1
}
agg_hostEP <- data.frame(unique(merged_dat$taxon_final),
                         agg_hostEP)

agg_hostEP[,2:length(agg_hostEP)] <- ifelse(agg_hostEP[,2:length(agg_hostEP)] > 0, 1, 0)
prevHost <- colSums(agg_hostEP[,2:length(agg_hostEP)])
summary(prevHost)
hist(prevHost)
dim(agg_hostEP)
table(prevHost > 50)

#DO FOR BACTERIA
library(tidyverse)
library(plyr)
library("wesanderson")


bacts <-  read.table("processedData/smallmem97_16S.sintax", fill = T)
bact_otus <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                     stringsAsFactors = F)
length(names(bact_otus)) -4

bact_otus[,3:length(bact_otus)] <- bact_otus[,3:length(bact_otus)] -1
bact_otus[,3:length(bact_otus)] <- bact_otus[,3:length(bact_otus)] / bact_otus$ISD

bact_meta <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                     stringsAsFactors = F)

bact_otus$sample <- gsub("X", "", bact_otus$sample)

merged_dat <- merge(bact_meta, bact_otus, by.x = "sample", by.y = "sample", all.y = T)
dim(merged_dat)
dim(bact_otus)
names(merged_dat)[200:205]

#recode the names of merged dat for each Zotu as the taxonomy hypothesis
for(i in 204:length(merged_dat)){
  if(any(grepl(names(merged_dat)[i], as.character(bacts$V1)))){
    names(merged_dat)[i] <- paste(names(merged_dat)[i],
                                  bacts$V4[bacts$V1 == names(merged_dat)[i]],
                                  sep = "_")
  }
}

#binarize
merged_dat[,grep("Zotu",names(merged_dat))] <- ifelse(merged_dat[,grep("Zotu",names(merged_dat))] > 0, 1, 0)

#agg by host 
hits <- list()
k <- 1
for(i in grep("Zotu",names(merged_dat))){
  hits[[k]] <- aggregate(merged_dat[,i] ~ merged_dat$taxon_final, FUN = sum)
  k <- k +1
}

#Overall prevalence
prevalence <- colSums(merged_dat[,grep("Zotu",names(merged_dat))], na.rm = T)
summary(prevalence)
prevalence[prevalence > 500]

#hist(prevalence)

#combine rows by host, rebinarize, and redetermine prevalence to see which things are in the most hosts
counts <- merged_dat[,grep("Zotu",names(merged_dat))]

agg_hostEN <- data.frame(matrix(nrow = length(unique(merged_dat$taxon_final)),
                                ncol = 1))
k <- 1
for(i in 1:length(counts)){
  agg_hostEN[,k] <- aggregate(counts[merged_dat$compartment == "EN",i] ~ merged_dat$taxon_final[merged_dat$compartment == "EN"], 
                              FUN = sum)[,2]
  colnames(agg_hostEN)[k] <- names(counts)[i]
  k <- k + 1
}
agg_hostEN <- data.frame(unique(merged_dat$taxon_final),
                         agg_hostEN)

agg_hostEN[,2:length(agg_hostEN)] <- ifelse(agg_hostEN[,2:length(agg_hostEN)] > 0, 1, 0)
prevHostENbact <- colSums(agg_hostEN[,2:length(agg_hostEN)])

#FOr EP
agg_hostEP <- data.frame(matrix(nrow = length(unique(merged_dat$taxon_final)),
                                ncol = 1))
k <- 1
for(i in 1:length(counts)){
  agg_hostEP[,k] <- aggregate(counts[merged_dat$compartment == "EP",i] ~ merged_dat$taxon_final[merged_dat$compartment == "EP"], FUN = sum)[,2]
  colnames(agg_hostEP)[k] <- names(counts)[i]
  k <- k + 1
}
agg_hostEP <- data.frame(unique(merged_dat$taxon_final),
                         agg_hostEP)

agg_hostEP[,2:length(agg_hostEP)] <- ifelse(agg_hostEP[,2:length(agg_hostEP)] > 0, 1, 0)
prevHostepbact <- colSums(agg_hostEP[,2:length(agg_hostEP)])
summary(prevHostepbact)


pdf(width = 8, height = 8, file = "visuals/occupancyCurves.pdf")
par(mfrow = c(2,2))
hist(prevHostEN, main = "fungal endophytes", las = 2, xlab = "Num. of Hosts")
hist(prevHost, main = "fungal epiphytes", las = 2, xlab = "Num. of Hosts")

hist(prevHostENbact, main = "bacterial endophytes", las = 2, xlab = "Num. of Hosts")

hist(prevHostepbact, main = "bacterial epiphytes", las = 2, xlab = "Num. of Hosts")
text(x = -15, y = 1600, xpd = NA, "d)", cex = 2)
text(x = -100, y = 1600, xpd = NA, "c)", cex = 2)
text(x = -15, y = 3900, xpd = NA, "b)", cex = 2)
text(x = -100, y = 3900, xpd = NA, "a)", cex = 2)
dev.off()

table(prevHostEN > length(unique(merged_dat$taxon_final))/2)
table(prevHost > length(unique(merged_dat$taxon_final))/2)
table(prevHostENbact > length(unique(merged_dat$taxon_final))/2)
table(prevHostepbact > length(unique(merged_dat$taxon_final))/2)


prevHostEN[prevHostEN == 62]
prevHost[prevHost == 62]
#cladosporidium was 140 and 3526 for fungi


#FIGURE OUT MOST PREVALENT MICROBE
rm(list=ls())
dat <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG_OCCUPANCY.csv")
head(dat)
dat <- dat[,-which(names(dat) %in% c("ISD", "duds"))]

observations <- NA
for(i in 3:length(dat)){
  observations[i] <- length(which(dat[,i] == 1))
}
tail(names(dat)[order(observations)])

table(dat$Zotu5710)
