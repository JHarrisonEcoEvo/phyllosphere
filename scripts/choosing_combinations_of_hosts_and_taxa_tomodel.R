#Choose combinations of host and microbial taxa to model based on microbial prevalence and how much I sampled
#the host taxon. 

rm(list=ls())

#load data
possibles <- read.csv("./processedData/ITS_taxa_to_model_via_randomforest.csv", stringsAsFactors = F)
X <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv", stringsAsFactors = F)
dat <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling", stringsAsFactors = F)

#Subset metadata to the most sampled hosts
tomodel <- names(which(table(X$taxon_final) > 30))
#Remove unknown stuff
tomodel <- tomodel[-grep("nknown", tomodel)]
X <- X[X$taxon_final %in% tomodel,]


##################################################
# wrangle the otu table as was done for modeling #
##################################################

#Note that we did considerable wrangling of these data using the
#"combine_pcr_dupes.*" script
dim(dat)
#remove taxa and samples with zero counts
dat <- dat[,colSums(dat[,3:length(dat)]) > 0]
dim(dat)
dat <- dat[rowSums(dat[,3:length(dat)]) > 0,]
dim(dat)

#remove rare microbes. THIS IS OK BC oF THE ISD, otherwise would be problematic
table(rowSums(dat[,3:length(dat)]) > 30)
dat <- dat[rowSums(dat[,3:length(dat)]) > 30,]

#extract the plant from the name and recode it so that each treatment
#is unique, use that to share info. 
treatments <- gsub("X(\\d+_\\d+_\\d+)_\\d+_(.*)","\\1\\2",names(dat))

labs <- seq(1,length(unique(treatments)),1)
k <- 1
for(i in unique(treatments)){
  treatments[treatments == i] <- paste(treatments[treatments == i], labs[k], sep = "_")
  k <- k + 1
}

#remove duplicates
dat <- dat[,-grep("dupe", treatments)]
treatments <- treatments[-grep("dupe", treatments)]

dim(dat)
dim(X)
dat[1:4,1:5]

#Need to get metadata into the same order as the modeled data
samps <- gsub("_16S_[1,2]_\\w+$", "", names(dat))
samps <- gsub("_16S$", "",samps)
samps <- gsub("dupe_", "",samps)
samps <- gsub("_rewash$", "",samps)
samps <- gsub("(E[NP]).*$", "\\1",samps)
samps <- gsub("(\\d+)(E[NP])$", "\\1_\\2",samps)
samps <- gsub("X", "",samps)

otus <- dat$otus
taxa <- data.frame(cbind(samps[3:length(samps)],
                   t(dat[,3:length(dat)])))
names(taxa) <- c("sampled",otus)
dim(taxa)
length(otus)

#OK, so the metadata should only have those host species that I sampled a lot. 

#Plan of attack for rest of script. 
# Merge the metadata and otu table
# Loop through by host species and find all OTUs that are fairly prevalent. 
# Write the host taxon and the OTU to two fields in an output file to feed into subsequent modeling. 

merged_dat <- merge(taxa, X, by.x = "sampled", by.y = "sample")
dim(merged_dat)
dim(taxa)
dim(X)

# Loop through by host species and find all OTUs that are fairly prevalent. 
keepers <- data.frame(matrix(ncol = 2, nrow = 1))
k <- 1
for(i in unique(merged_dat$taxon_final)){
  working <- merged_dat[merged_dat$taxon_final == i,2:(length(merged_dat)-198)]
  prevalence <- apply(working, 2, FUN=function(x){table(x>0)})
  for(j in 1:(length(prevalence)-2)){
    if(any(names(prevalence[[j]]) == "TRUE")){
      prev <- prevalence[[j]][which(names(prevalence[[j]])=="TRUE")]/sum(prevalence[[j]])
      if(prev >= 0.3){
        keepers[k, 1] <- names(prevalence)[j]
        keepers[k, 2] <- i
        k <- k + 1
      }
    }
  }  
}

write.csv(keepers, file = "./processedData/its_otus_to_analyze.csv")

##############
# DO FOR 16S #
##############

rm(list=ls())

#load data
possibles <- read.csv("./processedData/sixteenS_taxa_to_model_via_randomforest.csv", stringsAsFactors = F)
X <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv", stringsAsFactors = F)
dat <- read.csv("./processedData/otuTables/smallmem97_16s_for_modeling", stringsAsFactors = F)

#Subset metadata to the most sampled hosts
tomodel <- names(which(table(X$taxon_final) > 30))
#Remove unknown stuff
tomodel <- tomodel[-grep("nknown", tomodel)]
X <- X[X$taxon_final %in% tomodel,]

##################################################
# wrangle the otu table as was done for modeling #
##################################################

#Note that we did considerable wrangling of these data using the
#"combine_pcr_dupes.*" script
dim(dat)
#remove taxa and samples with zero counts
dat <- dat[,colSums(dat[,3:length(dat)]) > 0]
dim(dat)
dat <- dat[rowSums(dat[,3:length(dat)]) > 0,]
dim(dat)

#remove rare microbes. THIS IS OK BC oF THE ISD, otherwise would be problematic
table(rowSums(dat[,3:length(dat)]) > 30)
dat <- dat[rowSums(dat[,3:length(dat)]) > 30,]

#extract the plant from the name and recode it so that each treatment
#is unique, use that to share info. 
treatments <- gsub("X(\\d+_\\d+_\\d+)_\\d+_(.*)","\\1\\2",names(dat))

labs <- seq(1,length(unique(treatments)),1)
k <- 1
for(i in unique(treatments)){
  treatments[treatments == i] <- paste(treatments[treatments == i], labs[k], sep = "_")
  k <- k + 1
}

#remove duplicates
dat <- dat[,-grep("dupe", treatments)]
treatments <- treatments[-grep("dupe", treatments)]

dim(dat)
dim(X)
dat[1:4,1:5]

#Need to get metadata into the same order as the modeled data
samps <- gsub("_16S_[1,2]_\\w+$", "", names(dat))
samps <- gsub("_16S$", "",samps)
samps <- gsub("dupe_", "",samps)
samps <- gsub("_rewash$", "",samps)
samps <- gsub("(E[NP]).*$", "\\1",samps)
samps <- gsub("(\\d+)(E[NP])$", "\\1_\\2",samps)
samps <- gsub("X", "",samps)

otus <- dat$otus
taxa <- data.frame(cbind(samps[3:length(samps)],
                         t(dat[,3:length(dat)])))
names(taxa) <- c("sampled",otus)
dim(taxa)
length(otus)

#OK, so the metadata should only have those host species that I sampled a lot. 

#Plan of attack for rest of script. 
# Merge the metadata and otu table
# Loop through by host species and find all OTUs that are fairly prevalent. 
# Write the host taxon and the OTU to two fields in an output file to feed into subsequent modeling. 

merged_dat <- merge(taxa, X, by.x = "sampled", by.y = "sample")
dim(merged_dat)
dim(taxa)
dim(X)

# Loop through by host species and find all OTUs that are fairly prevalent. 
keepers <- data.frame(matrix(ncol = 2, nrow = 1))
k <- 1
for(i in unique(merged_dat$taxon_final)){
  working <- merged_dat[merged_dat$taxon_final == i,2:(length(merged_dat)-198)]
  prevalence <- apply(working, 2, FUN=function(x){table(x>0)})
  for(j in 1:(length(prevalence)-4)){
    if(any(names(prevalence[[j]]) == "TRUE")){
      prev <- prevalence[[j]][which(names(prevalence[[j]])=="TRUE")]/sum(prevalence[[j]])
      if(prev >= 0.3){
        keepers[k, 1] <- names(prevalence)[j]
        keepers[k, 2] <- i
        k <- k + 1
      }
    }
  }  
}

write.csv(keepers, file = "./processedData/sixteenS_otus_to_analyze.csv")




