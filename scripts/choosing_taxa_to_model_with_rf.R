rm(list=ls())
dat <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG", 
                fill = T, header = T, stringsAsFactors = F)

dat[,3:length(dat)] <- dat[,3:length(dat)] - 1

#convert data to qualitative data, 1 if present 0 otherwise
dat[,3:length(dat)][dat[,3:length(dat)] > 0] <- 1
dim(dat) #still good

#Figure out how many things were in 100 or more samples.
#This should be how many things I tried to model with the random forest
prev <- dat[,colSums(dat[,3:length(dat)]) >=100]

#Omit the isd, duds etc. 
prev <- prev[,grep("Zotu*", names(prev))]
dim(prev)

write.csv(names(prev), file = "./processedData/ITS_taxa_to_model_via_randomforest.csv", row.names = F)


############
# Bacteria
rm(list=ls())
dat <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG", 
                fill = T, header = T, stringsAsFactors = F)


dat[,3:length(dat)] <- dat[,3:length(dat)] - 1

#convert data to qualitative data, 1 if present 0 otherwise
dat[,3:length(dat)][dat[,3:length(dat)] > 0] <- 1
dim(dat) #still good

#Figure out how many things were in 100 or more samples.
#This should be how many things I tried to model with the random forest
prev <- dat[,colSums(dat[,3:length(dat)]) >=100]

#Omit the isd, duds etc. 
prev <- prev[,grep("Zotu*", names(prev))]
dim(prev)

write.csv(names(prev), file = "./processedData/sixteenS_taxa_to_model_via_randomforest.csv", row.names = F)
