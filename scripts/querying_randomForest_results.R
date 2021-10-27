rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/all_its_withhost.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

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
table(dat$mcc_nested_resampling)
table(dat$correctPositives)


dat <- read.csv("./modelingResults/all_its_occupancy.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$mcc_nested_resampling <- as.numeric(dat$mcc_nested_resampling)
summary(dat$mcc_nested_resampling)
table(dat$mcc_nested_resampling)
table(dat$correctPositives)


dat$taxon[which.max(dat$correctPositives)]
dat <- dat[
  order(dat$mcc_nested_resampling),]

tail(dat, n = 20)

#NO ISD TRANSFORMATION
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

