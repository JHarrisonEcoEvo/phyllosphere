rm(list=ls())

#Bring in results and figure out which taxa were predicted well enough to warrant
#extraction of important features. 

results <- read.csv("modelingResults/results_host//all_host_its_occupancy.csv", stringsAsFactors = F)
results <- results[results$taxon != "taxon",]
dim(results)
#results <- results[1:110,]
results <- results[as.numeric(results$mcc_nested_resampling) > 0.2,]
dim(results)

#Bring in variable importance metrics for only those taxa that were predicted
dat <- list.files("modelingResults/varImp_hostSpecific//")
dat <- dat[grep("*ITS*", dat)]
#dat <- dat[-grep("*rare*", dat)]
dat <- dat[grep("*OCC*", dat)]

dat <- grep(paste(results$taxon,"_", results$host, collapse="|", sep = ""), 
             dat, value=TRUE)

varimps <- lapply(paste("modelingResults/varImp_hostSpecific/", dat, sep = ""), read.csv)
length(varimps) == dim(results)[1]


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
 xtable(rev(sort(table(important)))[1:15])
 
 # % latex table generated in R 4.1.1 by xtable 1.8-4 package
 # % Fri Dec 10 15:51:06 2021
 # \begin{table}[ht]
 # \centering
 # \begin{tabular}{rr}
 # \hline
 # & important \\
 # \hline
 # SPAD\_420\_intensity &  42 \\
 # compartment\_EN &  42 \\
 # treeRich &  40 \\
 # Rel\_Chl\_intensity &  38 \\
 # gH. &  38 \\
 # phenology\_vegetative &  37 \\
 # densitometer &  37 \\
 # shrubRich &  36 \\
 # phenology\_fruiting &  33 \\
 # leaves\_extracted &  26 \\
 # slope\_perc &  16 \\
 # phenology\_flowering &   6 \\
 # julianDate &   6 \\
 # plant\_vol &   4 \\
 # deadDown &   4 \\
 # \hline
 # \end{tabular}
 # \end{table}
 # > 

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
