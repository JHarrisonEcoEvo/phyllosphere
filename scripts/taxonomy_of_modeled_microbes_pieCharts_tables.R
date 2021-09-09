#Make plots showing the taxonomy of the taxa that were modeled using random forest
# Groups of taxa to model include:
#1. 16s taxa modeled across all hosts
#2. ITS taxa modeled across all hosts
#1. 16s taxa modeled within hosts, by host. 
#2. ITS taxa modeled within hosts, by host.

#Extras. 
#which taxa were the most prevalent. Make a table. Were they modeled or not. 
#What was distribution of their abundances?

########################################
# 1. 16s taxa modeled across all hosts #
########################################

rm(list=ls())

tax <- read.csv("./processedData/sixteenS_taxa_to_model_via_randomforest.csv",stringsAsFactors = F)
head(tax)
sintax <- read.table("./processedData/smallmem97_16S.sintax", fill = T, stringsAsFactors = F)

hits <- sintax[sintax$V1 %in% tax$x,]
write.csv(hits, file = "./processedData/taxonomy_modeled_bacteria_landscape.csv", row.names = F)

#Extract unique phyla and classes
phyla <- gsub("k:Bacteria,p:(\\w+),.*","\\1",hits$V4)
phyla[phyla == "k:Bacteria"] <- "Unknown bacteria"
phyla[phyla == ""] <- "Unknown"
phyla <- gsub("k:Archaea,p:(\\w+),.*","\\1", phyla)
phyla[phyla == "k:Archaea"] <- "Unknown archaea"

class <- gsub("k:Bacteria,p:\\w+,c:(\\w+).*","\\1",hits$V4)
class <- gsub("k:Bacteria,p:\\w+,c:\\[(\\w+).*","\\1", class)
table(class)

#Pie chart
dat <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG")
dat[1:5,1:5] #1s were added!
dat[,3:length(dat)] <- dat[,3:length(dat)] - 1

#convert data to qualitative data, 1 if present 0 otherwise
dat[,3:length(dat)][dat[,3:length(dat)] > 0] <- 1
dim(dat) #still good

#Figure out how many things were in 100 or more samples.
#This should be how many things I tried to model with the random forest
prev <- dat[,colSums(dat[,3:length(dat)]) >=100]
prevalent <- dim(prev)[2] -2

hyperrare <- dim( dat[,colSums(dat[,3:length(dat)]) < 5])[2]
uncommon <- dim(dat[,colSums(dat[,3:length(dat)]) >= 5 & colSums(dat[,3:length(dat)]) < 100,])[2]

#how many were successfully modeled
results <- read.csv("processedData/all16s.csv")
results <- results[results$taxon != "taxon",]

good <- results[as.numeric(as.character(results$rsq_nested_resampling)) >= 0.01,]
predictable <- dim(good)[1]

#subtracting stuff so that the total of the pie chart is accurate, things not double counted
slice1 <- hyperrare
slice2 <- uncommon
slice3 <- prevalent - predictable
slice4 <- predictable

slices <- c(slice1, slice2, slice3, slice4)

#QC
sum(slices) == (dim(dat)-4)[2]
lbls <- c("Very rare OTUs  (< 5 samples)",
          "Uncommon OTUs (5-100 samples)",
          "Prevalent OTUs (100+ samples)",
          "Predictable OTUs (R2 > 1%)")

pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels 

pdf(width = 11, height = 11, file = "./visuals/16s_predictability_piechart.pdf")
pie(slices, labels = lbls, main="", col = c("gray","gray","gray","green"))
dev.off()

########################################
# 2. ITS taxa modeled across all hosts #
########################################

rm(list=ls())

tax <- read.csv("./processedData/ITS_taxa_to_model_via_randomforest.csv", stringsAsFactors = F)
head(tax)
dim(tax)
sintax <- read.table("./processedData/smallmem97_ITS.sintax", fill = T, stringsAsFactors = F)

hits <- sintax[sintax$V1 %in% tax$x,]
write.csv(hits, file = "./processedData/taxonomy_modeled_fungi_landscape.csv", row.names = F)

#Extract unique phyla and classes
phyla <- gsub("d:Fungi,p:(\\w+),.*","\\1",hits$V4)
phyla[phyla == "d:Fungi"] <- "Unknown bacteria"
phyla[phyla == ""] <- "Unknown"

class <- gsub("d:Fungi,p:\\w+,c:(\\w+).*","\\1",hits$V4)
class <- gsub("d:Fungi,p:\\w+,c:\\[(\\w+).*","\\1", class)
table(class)

#Pie chart
dat <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG")
dat[1:5,1:5] #1s were added!
dat[,3:length(dat)] <- dat[,3:length(dat)] - 1

#convert data to qualitative data, 1 if present 0 otherwise
dat[,3:length(dat)][dat[,3:length(dat)] > 0] <- 1
dim(dat) #still good

#Figure out how many things were in 100 or more samples.
#This should be how many things I tried to model with the random forest
prev <- dat[,colSums(dat[,3:length(dat)]) >=100]
prevalent <- dim(prev)[2] -2

hyperrare <- dim( dat[,colSums(dat[,3:length(dat)]) < 5])[2]
uncommon <- dim(dat[,colSums(dat[,3:length(dat)]) >= 5 & colSums(dat[,3:length(dat)]) < 100,])[2]

#how many were successfully modeled
results <- read.csv("processedData/allITS.csv")
results <- results[results$taxon != "taxon",]

good <- results[as.numeric(as.character(results$rsq_nested_resampling)) >= 0.01,]
predictable <- dim(good)[1]

#subtracting stuff so that the total of the pie chart is accurate, things not double counted
slice1 <- hyperrare
slice2 <- uncommon
slice3 <- prevalent - predictable
slice4 <- predictable

slices <- c(slice1, slice2, slice3, slice4)

#QC
sum(slices) == (dim(dat)-4)[2]
lbls <- c("Very rare OTUs  (< 5 samples)",
          "Uncommon OTUs (5-100 samples)",
          "Prevalent OTUs (100+ samples)",
          "Predictable OTUs (R2 > 1%)")

pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels 

pdf(width = 7, height = 7, file = "./visuals/its_predictability_piechart.pdf")
pie(slices, labels = lbls, main="", col = c("gray","gray","gray","green"))
dev.off()
########################################
# 1. 16s taxa modeled across all hosts #
########################################
