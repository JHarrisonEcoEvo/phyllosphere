rm(list=ls())

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

summary(as.numeric(dat$mse_nested_resampling))

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
