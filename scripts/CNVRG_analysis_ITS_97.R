#J. Harrison.
#Univ. of Wyoming.

rm(list=ls())
library(CNVRG)
#dat <- read.csv("./smallmem97_ITS_for_modeling",
#                fill = T, header = T, stringsAsFactors = F)
dat <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling",
                fill = T, header = T, stringsAsFactors = F)
tail(dat$otus)
#print dim to stdout
dim(dat)

#Note that we did considerable wrangling of these data using the
#"combine_pcr_dupes.*" script

#remove taxa and samples with zero counts
dat <- dat[,colSums(dat[,3:length(dat)]) > 0]
dim(dat)
dat <- dat[rowSums(dat[,3:length(dat)]) > 0,]
dim(dat)

#remove rare stuff
table(rowSums(dat[,3:length(dat)]) > 5)
dat <- dat[rowSums(dat[,3:length(dat)]) > 5,]
dim(dat)

#Remove the otus column and all columns before it
otus <- dat$otus 
dat <- dat[, -c(1,which(names(dat)=="otus"))]

#extract the plant from the name and recode it so that each treatment
#is unique, use that to share info. 
treatments <- gsub("X(\\d+_\\d+_\\d+)_\\d+_(.*)","\\1\\2",names(dat))

labs <- seq(1,length(unique(treatments)),1)
k <- 1
for(i in unique(treatments)){
  treatments[treatments == i] <- paste(treatments[treatments == i], labs[k], sep = "_")
  k <- k + 1
}

#sanity check
dim(dat)[2] == length(treatments)
cbind(names(dat), treatments)

#transpose the data and add a one
tdat <- t(dat)
tdat <- tdat+1
tdat <- data.frame(row.names(tdat), tdat)
names(tdat) <- c("sample", otus)
tdat$sample <- as.character(tdat$sample)

#sanity check
cbind(tdat$sample, treatments)
write.csv(tdat, file = "./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG")

write.csv(treatments, file = "./processedData/treatments_metadata.csv", row.names = F)

#convert ISD to integers. It is a float because of how I modified the ISD to reflect
#mistaken addition of too much ISD to certain plates. 
#See the combine_pcr_dupes... script
tdat$ISD <- round(tdat$ISD)

write.csv(unique(treatments), file = "./processedData/treatments_for_modeling_ITS.csv", 
          row.names = F)

modelOut <- cnvrg_VI(
  countData = tdat,
  starts = indexer(treatments)$starts,
  ends = indexer(treatments)$ends,
  #algorithm = "NUTS",
  #chains = 2,
  #burn = 500,
  #samples = 1500,
  #thinning_rate = 2,
  output_samples = 100,
  #  cores = 16,
  params_to_save = c("pi","p")
)
save.image(file = "/gscratch/jharri62/CNVRG_ITS_97vi.Rdata")

modelOut <- cnvrg_HMC(
  countData = tdat,
  starts = indexer(treatments)$starts,
  ends = indexer(treatments)$ends,
  algorithm = "NUTS",
  chains = 2,
  burn = 500,
  samples = 1500,
  thinning_rate = 2,
  #output_samples = 100,
    cores = 16,
  params_to_save = c("pi","p")
)
save.image(file = "/gscratch/jharri62/CNVRG_ITS_97vi_hmc.Rdata")

diffs <- diff_abund(model_output = modelOut, countData = tdat)
save(diffs, file = "ITS_diffs.Rdata")

ests <- extract_point_estimate(model_out = modelOut, countData = tdat)

write.csv(ests$pointEstimates_p , file = "multi_ests_its.csv")
write.csv(ests$pointEstimates_pi , file = "dir_ests_its.csv")

transformed <-
  isd_transform(model_out = modelOut, countData = tdat, isd_index = which(names(tdat) == "ISD"))

div <- diversity_calc(
  model_out = modelOut,
  countData = tdat[, 1:(length(tdat)-2)],
  params = "pi",
  entropy_measure = "shannon",
  equivalents = T
)
save(div, file = "ITS_div_isd_normalized.Rdata")
