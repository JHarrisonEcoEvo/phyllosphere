rm(list=ls())
dat <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
dat[,3:length(dat)] <- dat[,3:length(dat)] - 1
dat[,3:length(dat)] <- dat[,3:length(dat)] / dat$ISD

#names(dat)[length(dat)-4]
dat_occ <- dat[,3:(length(dat)-4)]
dat_occ <- data.frame(ifelse(dat_occ > 0, 1, 0))
#dat_occ[1:5,1:5]

cor.test(dat_occ[,1], dat[,3])

corout <- NA
for(i in 1:length(dat_occ)){
  corout[i] <- cor.test(dat_occ[,i], dat[,2+i])$estimate
}

summary(corout)

hist(corout)

#calculate median abundance versus prevalence (as proportion) for all taxa and correlate
#then do for just focal taxa

for(i in 3:length(dat)){
  dat[is.infinite(dat[,i]),i] <- NA
}

all_means <- colMeans(dat[,3:length(dat)], na.rm = T)

freq <- NA
for(i in 1:length(dat_occ)){
  if(any(names(table(dat_occ[,i] == 1)) == T)){
     tabout <- table(dat_occ[,i] == 1)
     freq[i] <- tabout[names(tabout) == "TRUE"] / sum(tabout)
  }else{
    freq[i] <- 0
  }
}

hist(freq)
summary(freq)

#length( all_means[length(all_means)-4])
all_means <- all_means[1:c(length(all_means)-4)]
cor.test(freq, 
         all_means)

# data:  freq and all_means
# t = 26.986, df = 2358, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.4543225 0.5160019
# sample estimates:
#   cor 
# 0.4857667 

focal <- read.csv("processedData/sixteenS_taxa_to_model_via_randomforest.csv")

cor.test(freq[names(all_means) %in% focal$x] , 
         all_means[names(all_means) %in% focal$x] )

# data:  freq[names(all_means) %in% focal$x] and all_means[names(all_means) %in% focal$x]
# t = 5.4164, df = 21, p-value = 2.258e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.5124941 0.8942257
# sample estimates:
#   cor 
# 0.7634224 

########Fungi
rm(list=ls())
dat <- read.csv("processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
dat[,3:length(dat)] <- dat[,3:length(dat)] - 1
dat[,3:length(dat)] <- dat[,3:length(dat)] / dat$ISD

#names(dat)[length(dat)-4]
dat_occ <- dat[,3:(length(dat)-4)]
dat_occ <- data.frame(ifelse(dat_occ > 0, 1, 0))
#dat_occ[1:5,1:5]

cor.test(dat_occ[,1], dat[,3])

corout <- NA
for(i in 1:length(dat_occ)){
  corout[i] <- cor.test(dat_occ[,i], dat[,2+i])$estimate
}

summary(corout)

hist(corout)

#calculate median abundance versus prevalence (as proportion) for all taxa and correlate
#then do for just focal taxa

for(i in 3:length(dat)){
  dat[is.infinite(dat[,i]),i] <- NA
}

all_means <- colMeans(dat[,3:length(dat)], na.rm = T)

freq <- NA
for(i in 1:length(dat_occ)){
  if(any(names(table(dat_occ[,i] == 1)) == T)){
    tabout <- table(dat_occ[,i] == 1)
    freq[i] <- tabout[names(tabout) == "TRUE"] / sum(tabout)
  }else{
    freq[i] <- 0
  }
}

hist(freq)
summary(freq)

#length( all_means[length(all_means)-4])
all_means <- all_means[1:c(length(all_means)-4)]
cor.test(freq, 
         all_means)
# data:  freq and all_means
# t = 31.365, df = 3185, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#   0.4588030 0.5118721
# sample estimates:
#   cor 
# 0.4857851 
focal <- read.csv("processedData/ITS_taxa_to_model_via_randomforest.csv")

cor.test(freq[names(all_means) %in% focal$x] , 
         all_means[names(all_means) %in% focal$x] )








