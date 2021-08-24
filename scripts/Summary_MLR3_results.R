
rm(list=ls())

#FIRST
#examine predictive ability of prevalent microbes at the landscape level. 

dat <- read.csv("./processedData/all16s.csv", stringsAsFactors = F)
head(dat)
#Remove the bogus rows
dat <- dat[-grep("taxon", dat$taxon),]

#Which can we predict. UGH.
#3 out of 50 
table(dat$rsq_nested_resampling > 0)
dat[dat$rsq_nested_resampling > 0,]

#Do for ITS
rm(list=ls())

dat <- read.csv("./processedData/allITS.csv", stringsAsFactors = F)
head(dat)
#Remove the bogus rows
dat <- dat[-grep("taxon", dat$taxon),]

#Which can we predict. 
#28/127
table(dat$rsq_nested_resampling > 0)
dat[dat$rsq_nested_resampling > 0,]

#SECOND
#examine predictive ability of prevalent microbes at the HOST level. 
#cant do anything. 

#begs question, is prediction at the landscape due to host or something else?

rm(list=ls())

dat <- read.csv("./processedData/all16s_host.csv", stringsAsFactors = F)
head(dat)
#Remove the bogus rows
dat <- dat[-grep("taxon", dat$taxon),]

table(dat$rsq_nested_resampling > 0)
dat[dat$rsq_nested_resampling > 0,]


rm(list=ls())
dat <- read.csv("./processedData/allITShost.csv", stringsAsFactors = F)
head(dat)
#Remove the bogus rows
dat <- dat[-grep("taxon", dat$taxon),]

table(dat$rsq_nested_resampling > 0)
dat[dat$rsq_nested_resampling > 0,]


######
#THIRD
# do with count data. to see if we get generally same results
#Same generally, but even worse.

#ITS is worse. 
rm(list=ls())
dat <- read.csv("./processedData/allITS_count.csv", stringsAsFactors = F)
head(dat)
#Remove the bogus rows
dat <- dat[-grep("taxon", dat$taxon),]

table(dat$rsq_nested_resampling > 0)
dat[dat$rsq_nested_resampling > 0,]

table(dat$rsq_reduced > 0)

rm(list=ls())
dat <- read.csv("./processedData/all16sCOUNT_landscape.csv", stringsAsFactors = F)
head(dat)
#Remove the bogus rows
dat <- dat[-grep("taxon", dat$taxon),]

table(dat$rsq_nested_resampling > 0)


