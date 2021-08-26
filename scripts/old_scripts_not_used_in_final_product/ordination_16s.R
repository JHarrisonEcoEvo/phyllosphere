rm(list = ls())
library(ecodist)
dat <- read.csv("./processedData/16sp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv",
                stringsAsFactors = F)
# dat <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
#                 stringsAsFactors = F)
#datr <-  vegan::drarefy(dat[,3:length(dat)], 2300)

datr <-  vegan::decostand(dat[,4:(length(dat)-4)], method = "hellinger")

datr[1:3,1:5]
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

# dat_bc <- bcdist(datr)
# ord <- pco(x = dat_bc)
#ord$vectors[,1] > 0.001
# plot(ord$vectors[,1], ord$vectors[,2], 
#      #type = "n", 
#      col = metadat$col,
#      xlab = "PCoA1", ylab = "PCoA2",
#      axes = TRUE, main = "PCoA (ecodist)")

dat_e <- dist(datr,method = "euclidean")
ord <- pco(x = dat_e)


#ord$vectors[,1] > 0.001
plot(ord$vectors[,1], ord$vectors[,2], 
     #type = "n", 
     col = metadat$col,
     xlab = "PCoA1", ylab = "PCoA2",
     axes = TRUE, main = "PCoA (ecodist)")
#what aer the weird ones 
metadat[which(ord$vectors[,1] < -0.05),]

#Try to figure out what the gradient is.
k <- 1
cor <- NA
name_cor <- NA
for(i in 1:length(metadat)){
        if(is.numeric(metadat[,i])){
                name_cor[k] <- names(metadat)[i]
               cor[k] <-   cor.test(ord$vectors[,2], metadat[,i])$estimate
               k <- k + 1
        }
}
# 
# ##################################
# # Do with most raw data possible #
# 
# rm(list = ls())
# library(ecodist)
# dat <- read.table("./processedData/otuTables/16s_97_clusterOTUs_otuTable",
#                 stringsAsFactors = F, header = T)
# 
# datr <-  vegan::drarefy(dat[,2:length(dat)], 2300)
# 
# datr[1:3,1:5]
# dim(dat)
# #This has ISD and all kinds of junk in it probably. Testing.
# dat_bc <- bcdist(datr)
#  
# ord <- pco(x = dat_bc)
# 
# pdf(width = 8, height = 8, file = "16sraw_ord.pdf")
# plot(ord$vectors[,1], ord$vectors[,2], 
#      xlab = "PCoA1", ylab = "PCoA2",
#      axes = TRUE, main = "PCoA (ecodist)")
# dev.off()



