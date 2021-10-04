rm(list=ls())

dat <- read.csv("./modelingResults/all16Slandscape.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

dat <- read.csv("./modelingResults/all16SlandscapeCount.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)




####
dat1 <- read.csv("./modelingResults/allITSlandscape.csv", stringsAsFactors = F)
head(dat1)
dat1 <- dat1[dat1$taxon != "taxon",]
dat1$rsq_nested_resampling <- as.numeric(dat1$rsq_nested_resampling)
summary(dat1$rsq_nested_resampling)

table(dat1$rsq_nested_resampling > 0.01)

dat <- read.csv("./modelingResults/allITSlandscapeCount.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

summary(as.numeric(dat$mse_nested_resampling))
summary(as.numeric(dat1$mse_nested_resampling))

dat1[dat1$rsq_nested_resampling > 0.01,]
dat[dat$rsq_nested_resampling > 0.01,]
summary(dat$rsq_nested_resampling )
