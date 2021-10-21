rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/all16S_landscape_count.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

dat <- read.csv("./modelingResults/allITS_landscape_count.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

dat$taxon[which.max(dat$rsq_nested_resampling)]

summary(as.numeric(dat$rsq_nested_resampling))

#host specific results
dat <- read.csv("./modelingResults/allITS_host_count.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)


dat <- read.csv("./modelingResults/all16s_host.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

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

tail(dat)
