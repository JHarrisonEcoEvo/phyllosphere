#Compare the models that used microbes as features to those that did not

############## ITS!
rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_landscape_isd/all_ITS.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dim(dat)
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)
dat[which.max(dat$rsq_nested_resampling),]
table(dat$rsq_nested_resampling > 0.25)


dat2 <- read.csv("./modelingResults/results_landscape_isd/allresultsITS_withmicrobes.csv", stringsAsFactors = F)
head(dat2)
dat2 <- dat2[dat2$taxon != "taxon",]
dim(dat2)
dat2$rsq_nested_resampling <- as.numeric(dat2$rsq_nested_resampling)
summary(dat2$rsq_nested_resampling)
length(dat2[,1])
table(dat2$rsq_nested_resampling > 0.01)
dat2[which.max(dat2$rsq_nested_resampling),]
table(dat2$rsq_nested_resampling > 0.25)


#Next directly compare each taxon

#Make sure the order is correct, it is. 
table(dat$taxon == dat2$taxon)

# Compare the R2 
# First set anything that is way low to zero. 
dat2$rsq_nested_resampling[dat2$rsq_nested_resampling < -1] <- 0
dat$rsq_nested_resampling[dat$rsq_nested_resampling < -1] <- 0

summary(dat2$rsq_nested_resampling - dat$rsq_nested_resampling)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.82534  0.01299  0.06531  0.07220  0.12061  0.78441 

##################
#Do for bacteria #
##################
rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_landscape_isd/all_16S.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dim(dat)
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)
dat[which.max(dat$rsq_nested_resampling),]
table(dat$rsq_nested_resampling > 0.25)


dat2 <- read.csv("./modelingResults/results_landscape_isd/allresults16S_withmicrobes.csv", stringsAsFactors = F)
head(dat2)
dat2 <- dat2[dat2$taxon != "taxon",]
dim(dat2)
dat2$rsq_nested_resampling <- as.numeric(dat2$rsq_nested_resampling)
summary(dat2$rsq_nested_resampling)
length(dat2[,1])
table(dat2$rsq_nested_resampling > 0.01)
dat2[which.max(dat2$rsq_nested_resampling),]
table(dat2$rsq_nested_resampling > 0.25)

#Next directly compare each taxon
#Make sure the order is correct, it is not.
table(dat$taxon == dat2$taxon)

#make it so
dat <- dat[dat$taxon %in% dat2$taxon,]
dat <- dat[!duplicated(dat$taxon),]

dat2 <- dat2[dat2$taxon %in% dat$taxon,]

dat2 <- dat2[!is.na(dat2$taxon),]

table(dat$taxon == dat2$taxon)


# Compare the R2 
# First set anything that is way low to zero. 
dat2$rsq_nested_resampling[dat2$rsq_nested_resampling < -1] <- 0
dat$rsq_nested_resampling[dat$rsq_nested_resampling < -1] <- 0

summary(dat2$rsq_nested_resampling - dat$rsq_nested_resampling)

#######################################
# Examine variable importance metrics #
# Curious to see how often other taxa are useful features #
#######################################
#16s
rm(list=ls())

#Bring in results and figure out which taxa were predicted well enough to warrant
#extraction of important features. 

results <- read.csv("./modelingResults/results_landscape_isd/allresults16S_withmicrobes.csv", stringsAsFactors = F)
results <- results[results$taxon != "taxon",]
results <- results[as.numeric(results$rsq_nested_resampling) > 0.01,]
dim(results)

#Bring in variable importance metrics for only those taxa that were predicted
dat <- list.files("modelingResults/varImp_16s_otherMicrobesAsFeatures/")

#find those files for the stuff that was predicted
dat <- grep(paste("variableImportance", results$taxon, "T16S.csv", sep = "",collapse="|"), 
            dat, value=TRUE)
#QC
length(dat) == dim(results)[1]

varimps <- lapply(paste("modelingResults/varImp_16s_otherMicrobesAsFeatures/", dat, sep = ""), read.csv)

#qc
length(dat) == length(varimps)# == dim(results)[1]

#Extract top ten from each and then do table to see which are most useful. 
#look at this clunker!

important <- vector()

for(i in 1:length(varimps)){
  important <- c(important,
                 as.character(varimps[[i]][head(rev(order(varimps[[i]][,2])), n = 10),1]))
}
table(important)
sort(table(important))
hist(sort(table(important)))
length(varimps)
#Wow other Zotus are way important.
#Extract the most important taxa
BestestBacteria <- names(sort(table(important)))[grep("Zotu",names(sort(table(important))))][9:18]

#Need to know their taxonomy and their average relative abundance
taxa <- read.csv("processedData/taxonomy_modeled_bacteria_landscape.csv")
output <- taxa[taxa$V1 %in% BestestBacteria, c(1,2)]

abunds <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG")
abunds[,3:length(abunds)] <- abunds[,3:length(abunds)] / abunds$ISD

abunds <- abunds[,names(abunds) %in% BestestBacteria]
abunds <- abunds[, match(output$V1,names(abunds))]

output$V1 == names(abunds)

output$median <- apply(abunds, 2, FUN=median)
output$mean <- apply(abunds, 2, FUN=mean)
output$max <- apply(abunds, 2, FUN=max)

xtable::xtable(output)

write.csv(output, file = "modelingResults/influential_bacterial_taxa.csv", row.names = F)

#####
#ITS#
#####

rm(list=ls())

#Bring in results and figure out which taxa were predicted well enough to warrant
#extraction of important features. 

results <- read.csv("./modelingResults/results_landscape_isd/allresultsITS_withmicrobes.csv", stringsAsFactors = F)
results <- results[results$taxon != "taxon",]
results <- results[as.numeric(results$rsq_nested_resampling) > 0.01,]
dim(results)

#Bring in variable importance metrics for only those taxa that were predicted
dat <- list.files("modelingResults/varImp_its_otherMicrobesAsFeatures/")

#find those files for the stuff that was predicted
dat <- grep(paste("variableImportance", results$taxon, "TITS.csv", sep = "",collapse="|"), 
            dat, value=TRUE)
#QC
length(dat) == dim(results)[1]

varimps <- lapply(paste("modelingResults/varImp_its_otherMicrobesAsFeatures/", dat, sep = ""), read.csv)

#qc
length(dat) == length(varimps)# == dim(results)[1]

#Extract top ten from each and then do table to see which are most useful. 
#look at this clunker!

important <- vector()

for(i in 1:length(varimps)){
  important <- c(important,
                 as.character(varimps[[i]][head(rev(order(varimps[[i]][,2])), n = 10),1]))
}
table(important)
sort(table(important))
hist(sort(table(important)))
length(varimps)
#Wow other Zotus are way important.

#Extract the most important taxa
BestestFungi <- names(sort(table(important)))[grep("Zotu",names(sort(table(important))))][148:157]

#Need to know their taxonomy and their average relative abundance
taxa <- read.csv("processedData/taxonomy_modeled_fungi_landscape.csv")
output <- taxa[taxa$V1 %in% BestestFungi, c(1,2)]

abunds <- read.csv("processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG")
abunds[,3:length(abunds)] <- abunds[,3:length(abunds)] / abunds$ISD

#See what the most abundant taxa are
overallMedians<- apply(abunds[,3:length(abunds)], 2, FUN=median)
head(rev(sort(overallMedians)))

abunds <- abunds[,names(abunds) %in% BestestFungi]
abunds <- abunds[, match(output$V1,names(abunds))]

output$V1 == names(abunds)

output$median <- apply(abunds, 2, FUN=median)
output$mean <- apply(abunds, 2, FUN=mean)
output$max <- apply(abunds, 2, FUN=max)

write.csv(output, file = "modelingResults/influential_fungal_taxa.csv", row.names = F)

xtable::xtable(output)
