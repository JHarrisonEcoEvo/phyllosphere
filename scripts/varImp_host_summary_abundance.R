rm(list=ls())

#Bring in results and figure out which taxa were predicted well enough to warrant
#extraction of important features. 

results <- read.csv("modelingResults/results_host//all_its_host.csv", stringsAsFactors = F)
dim(results)

results <- results[results$taxon != "taxon",]
dim(results)
results <- results[as.numeric(results$rsq_nested_resampling) > 0.01,]

#Bring in variable importance metrics for only those taxa that were predicted
dat <- list.files("modelingResults/varImp_hostSpecific////")
dat <- dat[grep("*ITS*", dat)]
dat <- dat[-grep("*educed*", dat)]
#dat <- dat[-grep("*NOHOST*", dat)]

dat <- grep(paste(results$taxon, collapse="|"), 
             dat, value=TRUE)

varimps <- lapply(paste("modelingResults/varImp_hostSpecific/", dat, sep = ""), read.csv)
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
xtable(rev(sort(table(important)))[1:15])
# % latex table generated in R 4.1.2 by xtable 1.8-4 package
# % Fri Nov 19 14:52:41 2021
# \begin{table}[ht]
# \centering
# \begin{tabular}{rr}
# \hline
# & important \\ 
# \hline
# SPAD\_420\_intensity &  25 \\ 
# densitometer &  25 \\ 
# compartment\_EN &  25 \\ 
# Rel\_Chl\_intensity &  24 \\ 
# phenology\_vegetative &  24 \\ 
# phenology\_fruiting &  22 \\ 
# gH. &  20 \\ 
# treeRich &  19 \\ 
# leaves\_extracted &  16 \\ 
# habit\_tree &  16 \\ 
# shrubRich &  11 \\ 
# slope\_perc &   8 \\ 
# habit\_shrub &   6 \\ 
# elev\_m &   5 \\ 
# deadDown &   5 \\ 
# \hline
# \end{tabular}
# \end{table}

# ###########################
# # Do for bacteria #####
#NO GOOD MODELS
# rm(list=ls())
# 
# #Bring in results and figure out which taxa were predicted well enough to warrant
# #extraction of important features. 
# 
#  results <- read.csv("modelingResults/results_host/all_16s_host_combinations_9.csv", stringsAsFactors = F)
# results <- results[results$taxon != "taxon",]
# 
# results <- results[results$rsq_nested_resampling > 0.01,]

# #Bring in variable importance metrics for only those taxa that were predicted
# dat <- list.files("modelingResults/varImp_hostCount/")
# dat <- dat[grep("*16S*", dat)]
# dat <- dat[-grep("*ITS*", dat)]
# 
# dat <- grep(paste(results$taxon, collapse="|"), 
#             dat, value=TRUE)
# 
# varimps <- lapply(paste("modelingResults/varImp_hostCount//", dat, sep = ""), read.csv)
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
