#############
# the fungi #
#############

rm(list=ls())
library(tidyverse)
library(plyr)
library("wesanderson")

its_otus <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                     stringsAsFactors = F)
length(names(its_otus)) -4

its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] -1
its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] / its_otus$ISD

its_meta <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                     stringsAsFactors = F)

its_otus$sample <- gsub("X", "", its_otus$sample)

merged_dat <- merge(its_meta, its_otus, by.x = "sample", by.y = "sample", all.y = T)

focal <- read.csv("processedData/ITS_taxa_to_model_via_randomforest.csv")

focal <- merged_dat[, names(merged_dat) %in% focal$x]

# var(aggregate(focal$Zotu103 ~ merged_dat$taxon_final, FUN = mean)[,2]) #order should be the same

# plot(type = "l",
#   aggregate(focal$Zotu103 ~ merged_dat$taxon_final, FUN = mean)[,2], ylim = c(0,10))
# 
# for(i in grep("Zotu", names(focal))){
#   lines(aggregate(focal[,i]~ merged_dat$taxon_final, FUN = mean)[,2])
# }

hits <- NA
k <-  1
for(i in grep("Zotu", names(focal))){
  hits[k]  <- length(unique(
   merged_dat$taxon_final[focal[,i] > 0]
 ))
 k <- k + 1
}

pdf(file = "visuals/hist_hostbreadth_focal.pdf", width = 12, height = 6)
par(mfrow = c(1,2))
hist(hits, main = "Freq. distribution of host bread for focal ITS", xlab = "Num of hosts")

#############
# the bacteria #
#############

#rm(list=ls())
library(tidyverse)
library(plyr)
library("wesanderson")

its_otus <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                     stringsAsFactors = F)
length(names(its_otus)) -4

its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] -1
its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] / its_otus$ISD

its_meta <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                     stringsAsFactors = F)

its_otus$sample <- gsub("X", "", its_otus$sample)

merged_dat <- merge(its_meta, its_otus, by.x = "sample", by.y = "sample", all.y = T)

focal <- read.csv("processedData/sixteenS_taxa_to_model_via_randomforest.csv")

focal <- merged_dat[, names(merged_dat) %in% focal$x]


# for(i in grep("Zotu", names(focal))){
#   lines(aggregate(focal[,i]~ merged_dat$taxon_final, FUN = mean)[,2])
# }

hits <- NA
k <-  1
for(i in grep("Zotu", names(focal))){
  hits[k]  <- length(unique(
    merged_dat$taxon_final[focal[,i] > 0]
  ))
  k <- k + 1
}

hist(hits, main = "Freq. distribution of host bread for focal 16S", xlab = "Num of hosts")
dev.off()
