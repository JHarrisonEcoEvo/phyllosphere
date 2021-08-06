rm(list = ls())
library(ecodist)
dat <- read.csv("./processedData/ITSp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv",
                stringsAsFactors = F)

metadat <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)
metadat$col <- "blue"
metadat$col[metadat$compartment == "EP"] <-"orange"

#match order
merged_dat <- merge(metadat, dat, by.x = "sample", by.y = "sample", all.y = T)
dim(merged_dat)

#################################################################
# Do ordinations for compartment and host taxon for all samples #
#################################################################

dat_bc <- bcdist(merged_dat[,204:(length(merged_dat)-2)])
#dat_e <- dist(normDat,method = "euclidean")
#hist(dat_e) #Way saturated
hist(dat_bc)
#ord <- pco(x = log(dat_e))
ord <- pco(x = dat_bc)

plot(ord$vectors[,1], ord$vectors[,2], 
     #type = "n", 
     col = merged_dat$col,
     pch = 16, 
     xlab = "PCoA1", ylab = "PCoA2",
     axes = TRUE, main = "PCoA (ecodist)")

merged_dat$col <- "blue"
for(i in unique(merged_dat$taxon.y)){
  merged_dat$col[merged_dat$taxon.y == i] <- randomcoloR::randomColor(count = 1)
}

plot(ord$vectors[,1], ord$vectors[,2], 
     #type = "n", 
     col = merged_dat$col,
     pch = 16, 
     xlab = "PCoA1", ylab = "PCoA2",
     axes = TRUE, main = "PCoA (ecodist)")


#####################################
# Perform ordinations by host taxon #
#####################################

merged_dat$col <- "blue"
merged_dat$col[merged_dat$compartment == "EP"] <-"orange"

dat_bc <- bcdist(merged_dat[merged_dat$taxon_final == "Artemisia tridentata" &
                              merged_dat$siteName == "pilot peak",
                            204:(length(merged_dat)-2)])

ord <- pco(x = dat_bc)
plot(ord$vectors[,1], ord$vectors[,2], 
     #type = "n", 
     col = merged_dat$col[merged_dat$taxon_final == "Artemisia tridentata"],
     pch = 16, 
     xlab = "PCoA1", ylab = "PCoA2",
     axes = TRUE, main = "PCoA (ecodist)")


############################
#sanity check doing with its90
############################
############################

dat <- read.table("./processedData/otuTables/its90_otuTable",
                  stringsAsFactors = F)
dat[1:3,1:5]

#Need to get metadata into the same order as the modeled data
sample <- gsub("_16S_[1,2]_\\w+$", "",dat[1,])
sample <- gsub("_16S$", "",sample)
sample <- gsub("dupe_", "",sample)
sample <- gsub("_rewash$", "",sample)
sample <- gsub("(E[NP]).*$", "\\1",sample)
sample <- gsub("(\\d+)(E[NP])$", "\\1_\\2",sample)

#Need to get metadata into the same order as the modeled data
sample <- gsub("_ITS_[1,2]_\\w+$", "",dat[1,])
sample <- gsub("ITS_[ATCG]*_[ATCG]*_", "",sample)
sample <- toupper(sample)
sample <- gsub("(\\d+)(E)", "\\1_\\2",sample)

metadat <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)

table(sample %in% metadat$sample) #some are missing, but good enough for qc
zotus <- dat[,1]

dat <- data.frame(t(dat[2:length(dat[,1]),2:length(dat)]))
dat[,1:3]
names(dat) <- zotus
dat$sample <- sample[2:length(sample)]

metadat$col <- "blue"
metadat$col[metadat$compartment == "EP"] <-"orange"

#match order
merged_dat <- merge(metadat, dat, by.x = "sample", by.y = "sample")
dim(merged_dat)
dim(dat)
dim(metadat)

#Note there are pcr duplicates. Just ignoring. 

dat_bc <- ecodist::bcdist(merged_dat[merged_dat$taxon_final == "Artemisia tridentata" &
                              merged_dat$siteName == "pilot peak",
                            204:(length(merged_dat)-2)])

ord <- ecodist::pco(x = dat_bc)
plot(ord$vectors[,1], ord$vectors[,2], 
     #type = "n", 
     col = merged_dat$col[merged_dat$taxon_final == "Artemisia tridentata"],
     pch = 16, 
     xlab = "PCoA1", ylab = "PCoA2",
     axes = TRUE, main = "PCoA (ecodist)")

#we get a similar pattern