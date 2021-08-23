rm(list = ls())
library(ecodist)
dat <- read.csv("./processedData/16sp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv",
                stringsAsFactors = F)
# dat <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
#                 stringsAsFactors = F)
#datr <-  vegan::drarefy(dat[,3:length(dat)], 2300)

dat[1:3,1:5]
tail(names(dat))
dim(dat)

metadat <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)
metadat$col <- "blue"
metadat$col[metadat$compartment == "EP"] <-"orange"

#sanity check
table(dat$sample == metadat$sample)

#sanity check that I use the correct last index during BC calculation
names(dat)[length(dat)-4]
# 
# #divide by ISD
# normDat <- dat[,5:(length(dat)-4)] / dat$ISD
# normDat[1:3,1:3]
#dat_bc <- bcdist(normDat)

dat_bc <- bcdist(dat[,4:(length(dat)-4)])
#dat_bc <- bcdist(datr[,3:2360])

#dat_e <- dist(normDat,method = "euclidean")

ord <- pco(x = dat_bc)

plot(ord$vectors[,1], ord$vectors[,2], 
     #type = "n", 
     col = metadat$col,
     xlab = "PCoA1", ylab = "PCoA2",
     axes = TRUE, main = "PCoA (ecodist)")
# 
# #Do over for stuff that is not super abundant as they are driving the ordination
# #which(colSums(normDat) > 100)
# normDat <- normDat[,colSums(normDat) < 100]
# dat_e <- dist(normDat,method = "euclidean")
# ord <- pco(x = log(dat_e))
# plot(ord$vectors[,1], ord$vectors[,2], 
#      #type = "n", 
#      col = metadat$col,
#      xlab = "PCoA1", ylab = "PCoA2",
#      axes = TRUE, main = "PCoA (ecodist)")
# 
# #Try looking at ratio of plant to bacteria
# hist(dat$plant)
# normDat <- dat[,5:(length(dat))] / dat$ISD
# normDat <- normDat / normDat$plant
# normDat[1:5,1:5]
# 
# dat_e <- bcdist(normDat[,1:(length(normDat)-4)])
# ord <- pco(x = dat_e)
# plot(ord$vectors[,1], ord$vectors[,2], 
#      #type = "n", 
#      col = metadat$col,
#      xlab = "PCoA1", ylab = "PCoA2",
#      axes = TRUE, main = "PCoA (ecodist)")
# 
# hist(rowSums(dat[,5:(length(dat)-4)]))
# 
# dat[dim(dat)[1],1:10]
