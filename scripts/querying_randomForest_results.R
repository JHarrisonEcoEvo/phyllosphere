
###############
#  landscape hellinger
##############

rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_landscape_hellinger/all_its_hella.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dim(dat)
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)
dat[which.max(dat$rsq_nested_resampling),]


dat <- read.csv("./modelingResults/results_landscape_hellinger/all_16s_hella.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)
 
###############
#  landscape isd
##############
rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_landscape_isd/all_ITS.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dim(dat)
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)
dat[which.max(dat$rsq_nested_resampling),]

rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_landscape_isd/all_ITS_NOHOST.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dim(dat)
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)
dat[which.max(dat$rsq_nested_resampling),]

# #Figure out median decline in R sq when comparing with and without host.
# dat <- read.csv("./modelingResults/results_landscape/all_ITS_withHostRetained_172", stringsAsFactors = F)
# dat <- dat[dat$taxon != "taxon",]
# dat$rsq_nested_resampling[as.numeric(dat$rsq_nested_resampling) < 0] <- 0
# 
# dat2 <- read.csv("./modelingResults/results_landscape/all_ITS_withHostremoved_172", stringsAsFactors = F)
# dat2 <- dat2[dat2$taxon != "taxon",]
# dat2$rsq_nested_resampling[as.numeric(dat2$rsq_nested_resampling) < 0] <- 0
# dim(dat)
# dim(dat2)
# 
# summary(as.numeric(dat$rsq_nested_resampling) - as.numeric(dat2$rsq_nested_resampling))
# 
# #make new df of the two different rsq and then subset to those rows that are greater than 0.01 for either
# #set. Then calculate decline summary stats
# rsqs <- data.frame(as.numeric(dat$rsq_nested_resampling),
#            as.numeric(dat2$rsq_nested_resampling))
# gooduns <- rsqs[rsqs[,1] > 0.01 | rsqs[,2] > 0.01,]
# summary(gooduns[,1] - gooduns[,2])

###############
# 16s landscape
##############


dat <- read.csv("./modelingResults/results_landscape_hellinger/all_16s_hella.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

hits <- dat$taxon

dat <- read.csv("./modelingResults/results_landscape_isd/all_16S.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)

rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_landscape_isd/all_16S_NOHOST.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dim(dat)
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.01)
dat[which.max(dat$rsq_nested_resampling),]


# dat <- read.csv("./modelingResults/results_landscape/all_16S_withHostremoved_23", stringsAsFactors = F)
# head(dat)
# dat <- dat[dat$taxon != "taxon",]
# dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
# summary(dat$rsq_nested_resampling)
# length(dat[,1])
# table(dat$rsq_nested_resampling > 0.01)
# 
# #Figure out median decline in R sq when comparing with and without host.
# dat <- read.csv("./modelingResults/results_landscape/all_16S_withHostRetained_23", stringsAsFactors = F)
# dat <- dat[dat$taxon != "taxon",]
# dat$rsq_nested_resampling[as.numeric(dat$rsq_nested_resampling) < 0] <- 0
# 
# dat2 <- read.csv("./modelingResults/results_landscape/all_16S_withHostremoved_23", stringsAsFactors = F)
# dat2 <- dat2[dat2$taxon != "taxon",]
# dat2$rsq_nested_resampling[as.numeric(dat2$rsq_nested_resampling) < 0] <- 0
# 
# summary(as.numeric(dat$rsq_nested_resampling) - as.numeric(dat2$rsq_nested_resampling))
# 
# #make new df of the two different rsq and then subset to those rows that are greater than 0.01 for either
# #set. Then calculate decline summary stats
# rsqs <- data.frame(as.numeric(dat$rsq_nested_resampling),
#                    as.numeric(dat2$rsq_nested_resampling))
# gooduns <- rsqs[rsqs[,1] > 0.01 | rsqs[,2] > 0.01,]
# summary(gooduns[,1] - gooduns[,2])

#######################
#landscape occupancy ITS
#####################
dat <- read.csv("./modelingResults/results_landscape_isd/all_ITS_OCC.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$mcc_nested_resampling <- as.numeric(dat$mcc_nested_resampling)
summary(dat$mcc_nested_resampling)
summary(as.numeric(dat$prop_positiveIDd))
table(as.numeric(dat$prop_positiveIDd) > 0.5)

table(dat$mcc_nested_resampling > 0.2)
#table(dat$mcc_nested_resampling)
#table(dat$correctPositives)

tail(dat[order(as.numeric(dat$prop_positiveIDd)),])

dat <- read.csv("./modelingResults/results_landscape_isd/all_16S_OCC.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dat$mcc_nested_resampling <- as.numeric(dat$mcc_nested_resampling)
summary(dat$mcc_nested_resampling)
table(dat$mcc_nested_resampling > 0.2)
#table(dat$correctPositives)
summary(as.numeric(dat$prop_positiveIDd))

summary(as.numeric(dat$correctPositives))

dat$taxon[which.max(dat$correctPositives)]
dat <- dat[
  order(dat$mcc_nested_resampling),]

tail(dat, n = 20)


##########################
# #NO ISD TRANSFORMATION, Hellinger standardized
##########################
# rm(list=ls())
# options(scipen = 99)
# dat <- read.csv("./modelingResults/results_landscapeNO_ISD/all_its_noisd", stringsAsFactors = F)
# head(dat)
# dat <- dat[dat$taxon != "taxon",]
# dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
# summary(dat$rsq_nested_resampling)
# length(dat[,1])
# table(dat$rsq_nested_resampling > 0.01)
# 
# rm(list=ls())
# options(scipen = 99)
# dat <- read.csv("./modelingResults/results_landscapeNO_ISD/all_16s_noisd", stringsAsFactors = F)
# head(dat)
# dat <- dat[dat$taxon != "taxon",]
# dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
# summary(dat$rsq_nested_resampling)
# length(dat[,1])
# table(dat$rsq_nested_resampling > 0.01)

#######
# ITS host combos
#########
rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_host/all_its_host.csv", stringsAsFactors = F)
dat <- dat[dat$taxon != "taxon",]
dim(dat)

table(dat$rsq_nested_resampling > 0.01)
#subset to the good ones
dat <- dat[dat$rsq_nested_resampling > 0.01,]
summary(as.numeric(dat$rsq_nested_resampling))

###
#Bacteria host combinations - all failed

rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_host/all_16s_host.csv", stringsAsFactors = F)
head(dat)
dim(dat)

dat <- dat[dat$taxon != "taxon",]
dim(dat)
table(dat$rsq_nested_resampling > 0.01)
#subset to the good ones
dat <- dat[dat$rsq_nested_resampling > 0.01,]
dat

# #Bacteria: No ISD, Helligner transformed, No host
# rm(list=ls())
# options(scipen = 99)
# dat <- read.csv("./modelingResults/results_landscapeNO_ISD//without_host/all_16s_nohost_noisd",
#                 stringsAsFactors = F)
# head(dat)
# 
# dat <- dat[dat$taxon != "taxon",]
# dim(dat)
# table(dat$rsq_nested_resampling > 0.01)
# 
# #Determine decline in performance. 
# dat_withHost <- read.csv("./modelingResults/results_landscapeNO_ISD/all_16s_noisd",
#                 stringsAsFactors = F)
# dat_withHost <- dat_withHost[dat_withHost$taxon != "taxon",]
# dat$taxon == dat_withHost$taxon
# 
# summary(as.numeric(dat$rsq_nested_resampling) - as.numeric(dat_withHost$rsq_nested_resampling))
# 
# #Fungi: No ISD, Helligner transformed, No host. FUNGI
# rm(list=ls())
# options(scipen = 99)
# dat <- read.csv("./modelingResults/results_landscapeNO_ISD/without_host/all_its_nohost_noisd",
#                 stringsAsFactors = F)
# head(dat)
# 
# dat <- dat[dat$taxon != "taxon",]
# dim(dat)
# table(dat$rsq_nested_resampling > 0.01)
# 
# #Determine decline in performance. 
# dat_withHost <- read.csv("./modelingResults/results_landscapeNO_ISD//all_its_noisd",
#                          stringsAsFactors = F)
# dat_withHost <- dat_withHost[dat_withHost$taxon != "taxon",]
# dat$taxon == dat_withHost$taxon
# 
# summary(as.numeric(dat$rsq_nested_resampling) - as.numeric(dat_withHost$rsq_nested_resampling))
# 
# dat[which.max(as.numeric(dat$rsq_nested_resampling) - as.numeric(dat_withHost$rsq_nested_resampling)),]

################################
#Host specific occupancy fungi
################################

rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_host/all_its_occupancy.csv",
                stringsAsFactors = F)
head(dat)

dat <- dat[dat$taxon != "taxon",]
dim(dat)
table(as.numeric(dat$mcc_nested_resampling) > 0.2)
summary(as.numeric(dat$mcc_nested_resampling ))
summary(as.numeric(dat$prop_positiveIDd ))

# #rarefied
# rm(list=ls())
# options(scipen = 99)
# dat <- read.csv("./modelingResults/results_host_occupancy/all_ITS_rare_110.csv",
#                 stringsAsFactors = F, fill = T)
# head(dat)
# 
# dat <- dat[dat$taxon != "taxon",]
# dim(dat)
# table(as.numeric(dat$mcc_nested_resampling) > 0.2)
# summary(as.numeric(dat$mcc_nested_resampling ))

################################
#Host specific occupancy bacteria
################################
rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_host/all_16s_occupancy.csv",
                stringsAsFactors = F)
head(dat)

dat <- dat[dat$taxon != "taxon",]
dim(dat)
table(as.numeric(dat$mcc_nested_resampling) > 0.2)
summary(as.numeric(dat$mcc_nested_resampling ))
summary(as.numeric(dat$prop_positiveIDd ))
# 
# #rarefied
# rm(list=ls())
# options(scipen = 99)
# dat <- read.csv("./modelingResults/results_host_occupancy/all_16S_rare_9.csv",
#                 stringsAsFactors = F, fill = T)
# head(dat)
# 
# dat <- dat[dat$taxon != "taxon",]
# dim(dat)
# table(as.numeric(dat$mcc_nested_resampling) > 0.2)
# summary(as.numeric(dat$mcc_nested_resampling ))