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


en <- merged_dat[merged_dat$compartment == "EN",
           grep("Zotu",names(merged_dat))]

en <- en[,colSums(en) > 0]
#ep
ep <- merged_dat[merged_dat$compartment == "EP",
                 grep("Zotu",names(merged_dat))]

ep <- ep[,colSums(ep) > 0]

colors <- wes_palette("FantasticFox1", n = 5, type = "discrete")

library(ggvenn)
pdf(width = 4, height = 4, file = "./visuals/fungi_venn.pdf")
ggvenn(list("EN" = names(en), "EP" = names(ep)), 
  stroke_size = 0.5, set_name_size = 5,
  fill_color  = c(colors[2], colors[3])
)
dev.off()

#############
# the bacteria #
#############

rm(list=ls())

bacts <- read.table("processedData/smallmem97_16S.sintax", fill = T)

bact_otus <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                      stringsAsFactors = F)
table(is.na(bact_otus[,3:length(bact_otus)] ))

bact_otus[,3:length(bact_otus)] <- bact_otus[,3:length(bact_otus)] - 1
bact_otus[,3:length(bact_otus)] <- bact_otus[,3:length(bact_otus)] / bact_otus$ISD
table(is.na(bact_otus[,3:length(bact_otus)] ))

bact_meta <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                      stringsAsFactors = F)

bact_otus$sample <- gsub("X", "", bact_otus$sample)

merged_dat <- merge(bact_meta, bact_otus, by.x = "sample", by.y = "sample", all.y = T)

#recode the names of merged dat for each Zotu as the taxonomy hypothesis
for(i in 204:length(merged_dat)){
  if(any(grepl(names(merged_dat)[i], bacts$V1))){
    names(merged_dat)[i] <- paste(names(merged_dat)[i],
                                  bacts$V4[bacts$V1 == names(merged_dat)[i]],
                                  sep = "_")
  }
}

merged_dat[1:10,1:10]

#Nas are from division by ISD I think
en <- merged_dat[merged_dat$compartment == "EN",
                 grep("Zotu",names(merged_dat))]

en <- en[,colSums(en, na.rm = T) > 0]
#ep
ep <- merged_dat[merged_dat$compartment == "EP",
                 grep("Zotu",names(merged_dat))]

ep <- ep[,colSums(ep, na.rm = T) > 0]
library(wesanderson)
colors <- wes_palette("FantasticFox1", n = 5, type = "discrete")

library(ggvenn)
pdf(width = 4, height = 4, file = "./visuals/bacteria_venn.pdf")
ggvenn(list("EN" = names(en), "EP" = names(ep)), 
       stroke_size = 0.5, set_name_size = 5,
       fill_color  = c(colors[1], colors[5])
)
dev.off()


