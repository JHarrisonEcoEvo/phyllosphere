rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/all_its_withhost.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)
dat[which.max(dat$rsq_nested_resampling),]
#Zotu6145


options(scipen = 99)
dat <- read.csv("./modelingResults/all_its_nohost.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

#Figure out median decline in R sq when comparing with and without host.
dat <- read.csv("./modelingResults/all_its_withhost.csv", stringsAsFactors = F)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling[as.numeric(dat$rsq_nested_resampling) < 0] <- 0

dat2 <- read.csv("./modelingResults/all_its_nohost.csv", stringsAsFactors = F)
dat2 <- dat2[dat2$taxon != "taxon",]
dat2$rsq_nested_resampling[as.numeric(dat2$rsq_nested_resampling) < 0] <- 0

summary(as.numeric(dat$rsq_nested_resampling) - as.numeric(dat2$rsq_nested_resampling))

#make new df of the two different rsq and then subset to those rows that are greater than 0.01 for either
#set. Then calculate decline summary stats
rsqs <- data.frame(as.numeric(dat$rsq_nested_resampling),
           as.numeric(dat2$rsq_nested_resampling))
gooduns <- rsqs[rsqs[,1] > 0.01 | rsqs[,2] > 0.01,]
summary(gooduns[,1] - gooduns[,2])

#plot(as.numeric(dat$rsq_nested_resampling) , as.numeric(dat2$rsq_nested_resampling))

#host specific results
dat <- read.csv("./modelingResults/allITS_host_count.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

########

dat <- read.csv("./modelingResults/all_16s_withhost.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

dat <- read.csv("./modelingResults/all_16s_nohost.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

#Figure out median decline in R sq when comparing with and without host.
dat <- read.csv("./modelingResults/all_16s_withhost.csv", stringsAsFactors = F)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling[as.numeric(dat$rsq_nested_resampling) < 0] <- 0

dat2 <- read.csv("./modelingResults/all_16s_nohost.csv", stringsAsFactors = F)
dat2 <- dat2[dat2$taxon != "taxon",]
dat2$rsq_nested_resampling[as.numeric(dat2$rsq_nested_resampling) < 0] <- 0

summary(as.numeric(dat$rsq_nested_resampling) - as.numeric(dat2$rsq_nested_resampling))

#make new df of the two different rsq and then subset to those rows that are greater than 0.01 for either
#set. Then calculate decline summary stats
rsqs <- data.frame(as.numeric(dat$rsq_nested_resampling),
                   as.numeric(dat2$rsq_nested_resampling))
gooduns <- rsqs[rsqs[,1] > 0.01 | rsqs[,2] > 0.01,]
summary(gooduns[,1] - gooduns[,2])


#occupancy
dat <- read.csv("./modelingResults/all_16s_occupancy.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$mcc_nested_resampling <- as.numeric(dat$mcc_nested_resampling)
summary(dat$mcc_nested_resampling)

table(dat$mcc_nested_resampling > 0.2)
table(dat$mcc_nested_resampling)
table(dat$correctPositives)


dat <- read.csv("./modelingResults/all_its_occupancy.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$mcc_nested_resampling <- as.numeric(dat$mcc_nested_resampling)
summary(dat$mcc_nested_resampling)
table(dat$mcc_nested_resampling > 0.2)
table(dat$correctPositives)


dat$taxon[which.max(dat$correctPositives)]
dat <- dat[
  order(dat$mcc_nested_resampling),]

tail(dat, n = 20)
# 
# #NO ISD TRANSFORMATION
rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_landscapeNO_ISD/all_its_noisd", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_landscapeNO_ISD/all_16s_noisd", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_host/all_its_byHostPlant.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
table(dat$rsq_nested_resampling > 0.01)
#subset to the good ones
dat <- dat[dat$rsq_nested_resampling > 0.01,]
dat
#check the sampling depth of host taxa
X <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv", stringsAsFactors = F)
X$samplename <- gsub("(\\d+_\\d+_\\d+_\\d+).*","\\1", X$samplename)
X <- X[!duplicated(X$samplename),]
sort(table(X$taxon_final))

dat$combo <- gsub(" ", "", 
                  paste(dat$taxon, dat$host, sep = "_"))
#bring in variable importance metrics

#Bring in variable importance metrics for only those taxa that were predicted
dat2 <- list.files("modelingResults/varImp_hostCount/")
dat2 <- dat2[grep("*ITS*", dat2)]
dat2 <- dat2[-grep("*Reduced*", dat2)]
#dat <- dat[-grep("*NOHOST*", dat)]

dat2 <- grep(paste(dat$combo, collapse="|"), 
            dat2, value=TRUE)

varimps <- lapply(paste("modelingResults/varImp_hostCount//", dat2, sep = ""), read.csv)
length(varimps)

#Extract top ten from each and then do table to see which are most useful. 

#look at this clunker!

important <- vector()

for(i in 1:length(varimps)){
  important <- c(important,
                 as.character(varimps[[i]][head(rev(order(varimps[[i]][,2])), n = 10),1]))
}
xtable::xtable(sort(table(important)),
               label = "table:hostFUN_important",
               caption = "Number of times a feature was in the top ten most important for models of fungal abundance that were limited to a single host taxon.")

sort(table(important))
hist(sort(table(important)))
length(varimps)




rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_host/all_16s_byHostPlant.csv", stringsAsFactors = F)
head(dat)

dat <- dat[dat$taxon != "taxon",]
dim(dat)
table(dat$rsq_nested_resampling > 0.01)
#subset to the good ones
dat <- dat[dat$rsq_nested_resampling > 0.01,]
dat
