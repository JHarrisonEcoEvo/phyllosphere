rm(list = ls())
library(ecodist)
dat <- read.csv("./processedData/ITSp_estimates_wrangled_for_post_modeling_analysis.csv",
                stringsAsFactors = F)
dat[1:3,1:5]
tail(names(dat))
dim(dat)

metadat <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)
metadat$col <- "blue"
metadat$col[metadat$compartment == "EP"] <-"orange"

#sanity check
table(dat$sample == metadat$sample)

#sanity check that I use the correct last index during BC calculation
#names(dat)[5]

#divide by ISD
normDat <- dat[,5:(length(dat)-2)] / dat$ISD
normDat[1:3,1:3]

dat_bc <- bcdist(normDat)
#dat_e <- dist(normDat,method = "euclidean")
#hist(dat_e) #Way saturated
hist(dat_bc)
#ord <- pco(x = log(dat_e))
ord <- pco(x = dat_bc)

plot(ord$vectors[,1], ord$vectors[,2], 
     #type = "n", 
     col = metadat$col,
     pch = 16, 
     xlab = "PCoA1", ylab = "PCoA2",
     axes = TRUE, main = "PCoA (ecodist)")

metadat$col <- "blue"
for(i in unique(metadat$taxon.y)){
  metadat$col[metadat$taxon.y == i] <- randomcoloR::randomColor(count = 1)
}

plot(ord$vectors[,1], ord$vectors[,2], 
     #type = "n", 
     col = metadat$col,
     pch = 16, 
     xlab = "PCoA1", ylab = "PCoA2",
     axes = TRUE, main = "PCoA (ecodist)")

