#J. Harrison.
#Univ. of Wyoming.

rm(list=ls())
library(CNVRG)
#dat <- read.csv("./smallmem97_ITS_for_modeling",
#                fill = T, header = T, stringsAsFactors = F)
dat <- read.csv("./processedData/smallmem97_ITS_for_modeling",
                fill = T, header = T, stringsAsFactors = F)

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
table(rowSums(dat[,3:length(dat)]) > 30)
dat <- dat[rowSums(dat[,3:length(dat)]) > 30,]

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
#remove duplicates
dat <- dat[,-grep("dupe", treatments)]
treatments <- treatments[-grep("dupe", treatments)]

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
write.csv(treatments, file = "./processedData/treatments_metadata.csv", row.names = F)

#convert ISD to integers. It is a float because of how I modified the ISD to reflect
#mistaken addition of too much ISD to certain plates. 
#See the combine_pcr_dupes... script
tdat$ISD <- round(tdat$ISD)

write.csv(unique(treatments), file = "./processedData/treatments_for_modeling_ITS.csv", 
          row.names = F)


#write.csv(tdat, row.names = F, file = "./processedData/ITS_otu_table_preModeling.csv")

table(is.numeric(tdat))
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
save.image(file = "/gscratch/jharri62/CNVRG_ITS_97.Rdata")


diffs <- diff_abund(model_output = modelOut, countData = tdat)
save(diffs, file = "ITS_diffs.Rdata")

ests <- extract_point_estimate(modelOut = modelOut, countData = tdat,treatments = length(unique(treatments)))
forExport <- data.frame(treatments, tdat[,1],ests$pointEstimates_p)
names(forExport)[2] <- "sample"
forExport <- data.frame(treatments, tdat[,1],ests$pointEstimates_pi)

write.csv(forExport, file = "ITS_pi_estimates.csv")

#Doing transformation...
load("CNVRG_ITS_97.Rdata")

transformed <-
  isd_transform(model_out = modelOut, countData = tdat, isd_index = 4267, format = "samples")

div <- diversity_calc(
  model_out = transformed$pi,
  countData = tdat[,(length(tdat)-2)],
  params = "pi",
  entropy_measure = "shannon",
  equivalents = T
)
save(div, file = "ITS_div.Rdata")


#sanity check, first part on teton
# vegdv <- vegan::diversity(tdat[,2:length(tdat)])
# boxplot(vegdv~treatments)
# treatments
# boxplot(vegdv~treatments)
# write.csv(data.frame(vegdv, treatments), file = "vegdv.csv")

#Looks fairly different without transforming/modeling
# test <- read.csv(file = "processedData/vegdv.csv")
# meta <- read.csv("processedData/treatments_metadata.csv")
# test$treatments <- gsub("(.*)_\\d+", "\\1", test$treatments)
# mtest <- merge(test, meta, by.x = "treatments", by.y = "treatmentClass")
# boxplot(mtest$vegdv ~ mtest$lifehistory + mtest$compartment)

# vegdv <- vegan::diversity(transformed) #gives same output as cnvrg, when using transformed data

rich0.001 <- rich_calc(model_out = modelOut, countData = tdat, params = "pi", threshold = 0.001)
rich0.005 <- rich_calc(model_out = modelOut, countData = tdat, params = "pi", threshold = 0.005)
rich0.01 <- rich_calc(model_out = modelOut, countData = tdat, params = "pi", threshold = 0.01)

