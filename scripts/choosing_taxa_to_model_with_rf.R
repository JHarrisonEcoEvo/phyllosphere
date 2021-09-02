rm(list=ls())
dat <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling", 
                fill = T, header = T, stringsAsFactors = F)

dat <- dat[rowSums(dat[,3:length(dat)]) > 0,]
total <- dim(dat)[1] #3908

#convert data to qualitative data, 1 if present 0 otherwise
dat[,3:length(dat)][dat[,3:length(dat)] > 0] <- 1
dim(dat) #still good

#Omit the isd, duds etc. 
dat <- dat[grep("Zotu*", dat$otus),]
dim(dat)

prevalent <- dat[rowSums(dat[,3:length(dat)]) > 100,]
dim(prevalent)
write.csv(prevalent$otus, file = "./processedData/ITS_taxa_to_model_via_randomforest.csv", row.names = F)


############
# Bacteria
rm(list=ls())
dat <- read.csv("./processedData/otuTables/smallmem97_16s_for_modeling", 
                fill = T, header = T, stringsAsFactors = F)

dat <- dat[rowSums(dat[,3:length(dat)]) > 0,]
total <- dim(dat)[1] 
dim(dat)

#convert data to qualitative data, 1 if present 0 otherwise
dat[,3:length(dat)][dat[,3:length(dat)] > 0] <- 1
dim(dat) #still good

#Omit the isd, duds etc. 
dat <- dat[grep("Zotu*", dat$otus),]
dim(dat)

prevalent <- dat[rowSums(dat[,3:length(dat)]) > 100,]
dim(prevalent)
write.csv(prevalent$otus, file = "./processedData/sixteenS_taxa_to_model_via_randomforest.csv", row.names = F)
