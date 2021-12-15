rm(list=ls())

#Bring in results and figure out which taxa were predicted well enough to warrant
#extraction of important features. 

results <- read.csv("modelingResults/results_host_hella/allits.csv", stringsAsFactors = F)
dim(results)

results <- results[results$taxon != "taxon",]
dim(results)
results <- results[as.numeric(results$rsq_nested_resampling) > 0.01,]
dim(results)
results$host <- gsub(" ","", results$host)

#Bring in variable importance metrics for only those taxa that were predicted
dat <- list.files("modelingResults/results_host_hella/")
dat <- dat[grep("HOSTvariableImportance*", dat)]
dat <- dat[grep("*ITS*", dat)]


dat <- grep(paste(results$taxon, "_",results$host, sep = "", collapse="|"), 
            dat, value=TRUE)

varimps <- lapply(paste("modelingResults/results_host_hella/", dat, sep = ""), read.csv)
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
