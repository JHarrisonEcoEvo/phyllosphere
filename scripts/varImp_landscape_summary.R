rm(list=ls())

#Bring in results and figure out which taxa were predicted well enough to warrant
#extraction of important features. 

results <- read.csv("modelingResults/allITS_landscape_count.csv", stringsAsFactors = F)
results <- results[results$rsq_nested_resampling > 0.01,]
results <- results[results$taxon != "taxon",]

#Bring in variable importance metrics for only those taxa that were predicted
dat <- list.files("modelingResults/varImp_landscapeCOUNT/")
dat <- dat[grep("*ITS*", dat)]
dat <- dat[-grep("*Reduced*", dat)]

dat <- grep(paste(results$taxon, collapse="|"), 
             dat, value=TRUE)

varimps <- lapply(paste("modelingResults/varImp_landscapeCOUNT/", dat, sep = ""), read.csv)
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

#Bring in taxonomy for each taxon and see if there are certain patterns by phyla
#where, say, a certain suite of predictors is more important for certain phyla
#If this appears to be the case, then make a heatmap figure


 




