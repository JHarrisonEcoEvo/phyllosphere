rm(list=ls())

#Bring in results and figure out which taxa were predicted well enough to warrant
#extraction of important features. 

results <- read.csv("modelingResults/results_host_occupancy/all_ITS_110.csv", stringsAsFactors = F)
results <- results[results$taxon != "taxon",]
dim(results)
results <- results[as.numeric(results$mcc_nested_resampling) > 0,]
dim(results)

#Bring in variable importance metrics for only those taxa that were predicted
dat <- list.files("modelingResults/varImp_hostCount_occupancy/")
dat <- dat[grep("*ITS*", dat)]
#dat <- dat[-grep("*rare*", dat)]
#dat <- dat[-grep("*Reduced*", dat)]

dat <- grep(paste(results$taxon, collapse="|"), 
             dat, value=TRUE)

varimps <- lapply(paste("modelingResults/varImp_hostCount_occupancy/", dat, sep = ""), read.csv)
length(varimps)

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

library(xtable)
xtable(rev(sort(table(important)))[1:10])

###########################
# Do for bacteria #####
#No good models
# rm(list=ls())
# 
# #Bring in results and figure out which taxa were predicted well enough to warrant
# #extraction of important features. 
# 
# results <- read.csv("modelingResults/results_host_occupancy/all_16S_9.csv", stringsAsFactors = F)
# results <- results[results$rsq_nested_resampling > 0.01,]
# results <- results[results$taxon != "taxon",]
# 
# #Bring in variable importance metrics for only those taxa that were predicted
# dat <- list.files("modelingResults/varImp_hostCount_occupancy//")
# dat <- dat[grep("*16S*", dat)]
# dat <- dat[-grep("*Reduced*", dat)]
# dat <- dat[-grep("*rare*", dat)]
# dat <- dat[grep("*CY.csv", dat)]
# 
# dat <- grep(paste(results$taxon, collapse="|"), 
#             dat, value=TRUE)
# 
# varimps <- lapply(paste("modelingResults/varImp_hostCount_occupancy//", dat, sep = ""), read.csv)
# length(varimps)
# 
# #Extract top ten from each and then do table to see which are most useful. 
# 
# #look at this clunker!
# 
# important <- vector()
# 
# for(i in 1:length(varimps)){
#   important <- c(important,
#                  as.character(varimps[[i]][head(rev(order(varimps[[i]][,2])), n = 10),1]))
# }
# table(important)
# sort(table(important))
# hist(sort(table(important)))
# length(varimps) #168
# 
# library(xtable)
# xtable(rev(sort(table(important)))[1:10])
