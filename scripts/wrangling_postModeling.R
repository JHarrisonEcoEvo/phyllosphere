#two part script one part for each locus
rm(list = ls())
dat <- read.csv("./processedData/16s_p_estimates.csv",
                stringsAsFactors = F)
dat[1:3,1:3]

metadat <- read.csv("./processedData/metadata_2.csv",
                    stringsAsFactors = F)
#CHANGE LOCUS AS NEEDED
metadat <- metadat[metadat$locus == "16S",]

dat$sample <- gsub("X", "",dat$sample)

#Need to get metadata into the same order as the modeled data
metadat$sample <- gsub("_ITS_[1,2]_\\w+$", "",metadat$sample)
metadat$sample <- gsub("_ITS$", "",metadat$sample)
metadat$sample <- gsub("dupe_", "",metadat$sample)
metadat$sample <- gsub("_rewash$", "",metadat$sample)
metadat$sample <- gsub("(E[NP]).*$", "\\1",metadat$sample)
metadat$sample <- gsub("(\\d+)(E[NP])$", "\\1_\\2",metadat$sample)

#Remove soil
metadat <- metadat[metadat$substrate == "plant",]

#There is no need to keep all PCR replicates from the metadata, since they have the 
#same metadata. Therefore, we will go through and keep just one row for each 
#unique sample identifier

metadat_reduced <- data.frame(matrix(nrow = length(unique(metadat$sample)),
                                     ncol = length(metadat)))
k <- 1
for(i in unique(metadat$sample)){
  metadat_reduced[k,] <- metadat[which(metadat$sample == i)[1],]
  k <- k + 1
}
names(metadat_reduced) <- names(metadat)


setdiff(metadat_reduced$sample, dat$sample)

#I spot checked some of these missing things and I removed them during wrangling becaaues
#They were cross-contaminated.

#make metadat and dat have same samples
metadat_reduced <- metadat_reduced[metadat_reduced$sample %in% dat$sample,]
#reorder so that the metadat matches the modeled dat
metadat_reduced <- metadat_reduced[match(dat$sample,
                                         metadat_reduced$sample),]

table(metadat_reduced$sample == dat$sample)

#Data are in proper order
#commence to doing analyses!
write.csv(metadat_reduced, file = "./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv")
write.csv(dat, file = "./processedData/16sp_estimates_wrangled_for_post_modeling_analysis.csv")

############
#Do for ITS#
############

rm(list = ls())
  dat <- read.csv("./processedData/ITS_multinomial_p_estimates.csv",
                stringsAsFactors = F)
dat[1:3,1:3]

metadat <- read.csv("./processedData/metadata_2.csv",
                    stringsAsFactors = F)
#CHANGE LOCUS AS NEEDED
metadat <- metadat[metadat$locus == "ITS",]

dat$sample <- gsub("X", "",dat$sample)

#Need to get metadata into the same order as the modeled data
metadat$sample <- gsub("_ITS_[1,2]_\\w+$", "",metadat$sample)
metadat$sample <- gsub("_ITS$", "",metadat$sample)
metadat$sample <- gsub("dupe_", "",metadat$sample)
metadat$sample <- gsub("_rewash$", "",metadat$sample)
metadat$sample <- gsub("(E[NP]).*$", "\\1",metadat$sample)
metadat$sample <- gsub("(\\d+)(E[NP])$", "\\1_\\2",metadat$sample)

#Remove soil
metadat <- metadat[metadat$substrate == "plant",]

#There is no need to keep all PCR replicates from the metadata, since they have the 
#same metadata. Therefore, we will go through and keep just one row for each 
#unique sample identifier

metadat_reduced <- data.frame(matrix(nrow = length(unique(metadat$sample)),
                                     ncol = length(metadat)))
k <- 1
for(i in unique(metadat$sample)){
  metadat_reduced[k,] <- metadat[which(metadat$sample == i)[1],]
  k <- k + 1
}
names(metadat_reduced) <- names(metadat)


setdiff(metadat_reduced$sample, dat$sample)

#I spot checked some of these missing things and I removed them during wrangling becaaues
#They were cross-contaminated.

#make metadat and dat have same samples
metadat_reduced <- metadat_reduced[metadat_reduced$sample %in% dat$sample,]
#reorder so that the metadat matches the modeled dat
metadat_reduced <- metadat_reduced[match(dat$sample,
                                         metadat_reduced$sample),]

table(metadat_reduced$sample == dat$sample)

metadat_reduced$samplename[is.na(metadat_reduced$compartment)] 

metadat_reduced$compartment[is.na(metadat_reduced$compartment)] <- "EN"

#Data are in proper order
#commence to doing analyses!
write.csv(metadat_reduced, file = "./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv")
write.csv(dat, file = "./processedData/ITSp_estimates_wrangled_for_post_modeling_analysis.csv")
