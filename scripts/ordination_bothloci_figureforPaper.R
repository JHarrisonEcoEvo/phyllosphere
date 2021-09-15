rm(list = ls())
library(ecodist)
library(vegan)
library(dendextend)

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#################################################################
# Do ordinations for compartment and host taxon for all samples #
#################################################################

#################################################################
# ITS processing  #
#################################################################

dat <- read.csv("./processedData/ITSp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv",
                stringsAsFactors = F)

metadat <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)

#match order
merged_dat <- merge(metadat, dat, by.x = "sample", by.y = "sample", all.y = T)
dim(merged_dat)
dim(dat)
dim(metadat)
tail(names(merged_dat))
head(names(merged_dat),205)
unique(merged_dat$habit)

datr <-  vegan::decostand(merged_dat[,202:(length(merged_dat)-2)], method = "hellinger")

dat_e <- dist(datr,method = "euclidean")

ord <- cmdscale(dat_e, eig = TRUE, k = 2)

## Assign colors for ordination
library(wesanderson)
colors <- wes_palette("FantasticFox1", n =  4, type = "discrete")

merged_dat$col <- colors[1]
merged_dat$col[merged_dat$compartment == "EP"] <- colors[3]

## See if the ordination correlates with total abundance or the abundance of the duds
## the latter where the unidentifiable sequences
cor.test(ord$points[,2], rowSums(merged_dat[,202:(length(merged_dat)-2)]))
cor.test(ord$points[,2], merged_dat$duds)

###################
# ITS EN vs EP ASVs  #
###################

pdf(width = 8, height = 8, file = "./visuals/ordination_ITS_enVsep_asvs_multi_div_by_ISD.pdf")

plot(ord$points[,1], ord$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat$col, 0.8),
     pch = 16,
     main = "",
     ylim = c(-0.7, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
text(x = -0.8, y=0.55, "a)", xpd = NA, cex = 1.7)

#ordihull(ord, merged_dat$col, col = unique(merged_dat$col))
legend(x = -0.2, y = 0.5, 
       legend = c("EN","EP"), 
       col = unique(add.alpha(merged_dat$col, 0.8)), 
       pch = 19, 
       horiz = T,
       cex = 2, 
       xpd = NA,
       bty = "n")

dev.off()
######################################

## Redo the colors so that they are by host taxon
#colors <- randomcoloR::randomColor(count = length(unique(merged_dat$taxon_final)))
colors <- rainbow(n = length(unique(merged_dat$taxon_final)))
k <- 1
for(i in unique(merged_dat$taxon_final)){
  merged_dat$col[merged_dat$taxon_final == i] <- colors[k]
  k <- k + 1
}

###################
# ITS by host taxon - ASVs  #
###################

pdf(width = 8, height = 8, file = "./visuals/ordination_ITS_byHost_asvs_multi_div_by_ISD.pdf")

plot(ord$points[,1], ord$points[,2],
     #type = "n",
     col = add.alpha(merged_dat$col, 0.8),
     pch = 16,
     main = "",
     ylim = c(-0.7, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
text(x = -0.8, y=0.55, "a)", xpd = NA, cex = 1.7)

#Ordihull is not working correctly. Do not feel like debugging.
#ordihull(ord, merged_dat$taxon_final, col = unique(merged_dat$col))
dev.off()


# plot(ord$points[,1], ord$points[,2],
#      #type = "n",
#      col = add.alpha(merged_dat$col, 0.8),
#      pch = 16,
#      main = "",
#      xlim = c(-0.05,0.1),
#      ylim = c(0, 0.1),
#      cex = 1.3,
#      frame.plot = F,
#      xlab = "PCoA1", ylab = "PCoA2",
#      yaxt = "n")
# axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
# text(x = -0.1, y=0.6, "b)", xpd = NA, cex = 1.7)

###################
# Legend plot #
###################

pdf(width = 7, height = 14, file = "./visuals/ordination_hostLEGEND.pdf")

plot(ord$points[,1], ord$points[,2], type = "n",axes = F, ylab = "", xlab = "")
legend(x = 0, y = 0.2,
       legend = unique(merged_dat$taxon_final),
       col = unique(add.alpha(merged_dat$col, 0.8)),
       pch = 19,
       cex = 0.6,
       xpd = NA,
       bty = "n")

dev.off()

###################
# ITS by life history  #
###################

#colors <- rainbow(n = length(unique(merged_dat$lifehistory)))
library(wesanderson)
colors <- wes_palette("FantasticFox1", n =  8, type = "continuous")
colors <- colors[c(1,3,4,7)]

k <- 1
for(i in unique(merged_dat$habit)){
  merged_dat$col[merged_dat$habit == i] <- colors[k]
  k <- k + 1
}

pdf(width = 8, height = 8, file = "./visuals/ordination_ITS_lifehistory_asvs_multi_div_by_ISD.pdf")

plot(ord$points[,1], ord$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat$col, 0.8),
     pch = 16,
     main = "",
     ylim = c(-0.7, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
text(x = -0.8, y=0.55, "b)", xpd = NA, cex = 1.7)

legend(x = -.10, y = -0.4, 
       legend = unique(merged_dat$habit), 
       col = unique(add.alpha(merged_dat$col, 0.8)), 
       pch = 19, 
       cex = 1.8, 
       xpd = NA,
       bty = "n")
#Ordihull is not working correctly. Do not feel like debugging.
#ordihull(ord, merged_dat$taxon_final, col = unique(merged_dat$col))
dev.off()

###########################
# bacteria data wrangling #
###########################
rm(list = ls())
library(ecodist)
library(vegan)
library(dendextend)

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

dat16s <- read.csv("./processedData/16sp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv",
                stringsAsFactors = F)

metadat16s <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)

#match order
merged_dat16s <- merge(metadat16s, dat16s, by.x = "sample", by.y = "sample", all.y = T)
dim(merged_dat16s)

head(names(merged_dat16s),206)
tail(names(merged_dat16s),206)
names(merged_dat16s)[length(merged_dat16s)-4]

datr16s <-  vegan::decostand(merged_dat16s[,202:(length(merged_dat16s)-4)], method = "hellinger")

dat_e16s <- dist(datr16s, method = "euclidean")

ord16s <- cmdscale(dat_e16s, eig = TRUE, k = 2)

library(wesanderson)
colors <- wes_palette("FantasticFox1", n =  4, type = "discrete")

merged_dat16s$col <- colors[1]
merged_dat16s$col[merged_dat16s$compartment == "EP"] <- colors[3]
###########################
# bacteria EP vs EN ordination #
###########################

pdf(width = 8, height = 8, file = "./visuals/ordination_16S_epVsen_asvs_hellingerEuclidean_multi_div_by_ISD.pdf")

plot(ord16s$points[,1], ord16s$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat16s$col, 0.8),
     pch = 16,
     main = "",
     xlim = c(-0.6, 0),
     ylim = c(-1, 0.8),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.8,0.4,0.4), labels =  seq(-0.8,0.4,0.4), las = 2)
text(x = -0.6, y=0.6, "e)", xpd = NA, cex = 1.7)

#ordihull(ord, merged_dat$col, col = unique(merged_dat$col))
legend("topright",
       legend = c("EN","EP"), 
       col = unique(add.alpha(merged_dat16s$col, 0.8)), 
       pch = 19, 
       horiz = T,
       cex = 2, 
       xpd = NA,
       bty = "n")

dev.off()


###########################
# bacteria EP vs EN ordination zoomed in #
###########################

pdf(width = 8, height = 8, file = "./visuals/ordination_16S_epVsen_asvs_zoomedInhellingerEuclidean_multi.pdf")

plot(ord16s$points[ord16s$points[,1] > -0.02,1], ord16s$points[ord16s$points[,1] > -0.02,2], 
     #type = "n", 
     col = add.alpha(merged_dat16s$col, 0.5),
     pch = 16,
     main = "",
     ylim = c(-0.1, 0.1),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.1,0.1,0.05), labels =   seq(-0.1,0.1,0.05), xpd = NA, las = 2)
#ordihull(ord, merged_dat$col, col = unique(merged_dat$col))
legend(x = -0.01, y = 0.1, 
       legend = c("EN","EP"), 
       col = unique(add.alpha(merged_dat16s$col, 0.5)), 
       pch = 19, 
       horiz = T,
       cex = 2, 
       xpd = NA,
       bty = "n")

dev.off()

#########################
# Redoing the colors
# 
# colors <- rainbow(n = length(unique(merged_dat16s$taxon_final)))
# k <- 1
# for(i in unique(merged_dat16s$taxon_final)){
#   merged_dat16s$col[merged_dat16s$taxon_final == i] <- colors[k]
#   k <- k + 1
# }
# 
# ###########################
# # bacteria by host - asvs #
# ###########################
# 
# pdf(width = 8, height = 8, file = "./visuals/ordination_16S_host_asvshellingerEuclidean_multi.pdf")
# 
# plot(ord16s$points[,1], ord16s$points[,2], 
#      #type = "n", 
#      col = add.alpha(merged_dat16s$col, 0.5),
#      pch = 16,
#      main = "",
#      ylim = c(-0.7, 0.4),
#      cex = 1.3,
#      frame.plot = F,
#      xlab = "PCoA1", ylab = "PCoA2",
#      yaxt = "n")
# axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
# text(x = -0.1, y=0.6, "f)", xpd = NA, cex = 1.7)
# 
# dev.off()
# 
# #Zoomed in
# pdf(width = 8, height = 8, file = "./visuals/ordination_16S_host_asvsZoomedinhellingerEuclidean_multi.pdf")
# 
# plot(ord16s$points[ord16s$points[,1] > -0.04,1], ord16s$points[ord16s$points[,1] > -0.04,2], 
#      #type = "n", 
#      col = add.alpha(merged_dat16s$col, 0.5),
#      pch = 16,
#      main = "",
#      ylim = c(-0.7, 0.4),
#      cex = 1.3,
#      frame.plot = F,
#      xlab = "PCoA1", ylab = "PCoA2",
#      yaxt = "n")
# axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
# text(x = -0.1, y=0.6, "f)", xpd = NA, cex = 1.7)
# 
# dev.off()


###################
# 16 by life history  #
###################

#colors <- rainbow(n = length(unique(merged_dat$lifehistory)))
library(wesanderson)
colors <- wes_palette("FantasticFox1", n =  8, type = "continuous")
colors <- colors[c(1,3,4,7)]

k <- 1
for(i in unique(merged_dat16s$habit)){
  merged_dat16s$col[merged_dat16s$habit == i] <- colors[k]
  k <- k + 1
}

pdf(width = 8, height = 8, file = "./visuals/ordination_16S_lifehistory_asvshellingerEuclidean_multi.pdf")

plot(ord16s$points[,1], ord16s$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat16s$col, 0.8),
     pch = 16,
     main = "",
     ylim = c(-0.7, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
options(scipen = 99)
axis(side = 2, at = c(-0.6,-0.4,-.2, 0, .2, .4), labels =   c(-0.6,-0.4,-.2, 0, .2, .4), las = 2)
text(x = -0.1, y=0.6, "b)", xpd = NA, cex = 1.7)

legend(x = -0.1, y = -0.4, 
       legend = unique(merged_dat16s$habit), 
       col = unique(add.alpha(merged_dat16s$col, 0.8)), 
       pch = 19, 
       cex = 2, 
       xpd = NA,
       bty = "n")

dev.off()


pdf(width = 8, height = 8, file = "./visuals/ordination_16S_lifehistory_asvsZoomedhellingerEuclidean_multi.pdf")

plot(ord16s$points[,1], ord16s$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat16s$col, 0.8),
     pch = 16,
     main = "",
     xlim = c(-0.02, 0.02),
     ylim = c(-0.02, 0.04),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")

axis(side = 2, at = c(-0.02,-.01, 0, .01, .02,.03,.04), labels =c(-0.02,-.01, 0, .01, .02,.03,.04), las = 2)
text(x = -0.1, y=0.6, "b)", xpd = NA, cex = 1.7)

legend(x = 0.01, y = 0.04, 
       legend = unique(merged_dat16s$habit), 
       col = unique(add.alpha(merged_dat16s$col, 0.8)), 
       pch = 19, 
       cex = 2, 
       xpd = NA,
       bty = "n")

dev.off()

############################
# ordination of count data #
############################
#TRYING UNNORMALIZED DATA

rm(list = ls())
library(ecodist)
library(vegan)
library(dendextend)

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#################################################################
# Do ordinations for compartment and host taxon for all samples #
#################################################################

#################################################################
# ITS processing  #
#################################################################

dat <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
dat$sample <- gsub("X","", dat$sample)
metadat <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)

#match order
merged_dat <- merge(metadat, dat, by.x = "sample", by.y = "sample", all.y = T)
dim(merged_dat)
dim(dat)
dim(metadat)
tail(names(merged_dat))
head(names(merged_dat),205)
unique(merged_dat$lifehistory)

#datr <-  vegan::decostand(merged_dat[,202:(length(merged_dat)-2)], method = "hellinger")
datr <-  vegan::decostand(merged_dat[,202:(length(merged_dat)-2)], method = "log")

dat_e <- vegan::vegdist(datr, method =  "altGower")

ord <- cmdscale(dat_e, eig = TRUE, k = 2)

## Assign colors for ordination
library(wesanderson)
colors <- wes_palette("FantasticFox1", n =  4, type = "discrete")

merged_dat$col <- colors[1]
merged_dat$col[merged_dat$compartment == "EP"] <- colors[3]

## See if the ordination correlates with total abundance or the abundance of the duds
## the latter where the unidentifiable sequences
cor.test(ord$points[,2], rowSums(merged_dat[,202:(length(merged_dat)-2)]))
cor.test(ord$points[,2], merged_dat$duds)

###################
# ITS EN vs EP ASVs  #
###################
pdf(width = 8, height = 8, file = "./visuals/ordination_ITS_enVsep_asvsGowerCount_no_isd.pdf")

plot(ord$points[,1], ord$points[,2],
     #type = "n",
     col = add.alpha(merged_dat$col, 0.8),
     pch = 16,
     main = "",
     ylim = c(-0.4, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.4,0.4,0.2), labels =  seq(-0.4,0.4,0.2), las = 2)
text(x = -0.4, y=0.4, "a)", xpd = NA, cex = 1.7)

#ordihull(ord, merged_dat$col, col = unique(merged_dat$col))
legend(x = 0.01, y = 0.4, 
       legend = c("EN","EP"),
       col = add.alpha(c(colors[1],colors[3]), 0.8), 
       pch = 19, 
       cex = 2, 
       xpd = NA,
       bty = "n")
dev.off()
######################################

###################
# ITS by life history  #
###################

#colors <- rainbow(n = length(unique(merged_dat$lifehistory)))
library(wesanderson)
colors <- wes_palette("FantasticFox1", n =  8, type = "continuous")
colors <- colors[c(1,3,4,7)]

k <- 1
for(i in unique(merged_dat$lifehistory)){
  merged_dat$col[merged_dat$lifehistory == i] <- colors[k]
  k <- k + 1
}

pdf(width = 8, height = 8, file = "./visuals/ordination_ITS_lifehistory_asvs_Gower_no_isd.pdf")

plot(ord$points[,1], ord$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat$col, 0.8),
     pch = 16,
     main = "",
     ylim = c(-0.4, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.4,0.4,0.2), labels =  seq(-0.4,0.4,0.2), las = 2)
text(x = -0.4, y=0.4, "a)", xpd = NA, cex = 1.7)

legend(x = 0.01, y = 0.5, 
       legend = unique(merged_dat$habit),
       col = unique(add.alpha(merged_dat$col, 0.8)),
       pch = 19,
       #horiz = T,
       cex = 2,
       xpd = NA,
       bty = "n")
dev.off()


####Divide by ISD 
merged_dat[,202:(length(merged_dat)-2)] <- merged_dat[,202:(length(merged_dat)-2)]/merged_dat$ISD

datr <-  vegan::decostand(merged_dat[,202:(length(merged_dat)-2)], 
                          method = "hellinger")

dat_e <- vegan::vegdist(datr, method =  "euclidean")

ord <- cmdscale(dat_e, eig = TRUE, k = 2)

pdf(width = 8, height = 8, file = "./visuals/ordination_ITS_lifehistory_isdCounts_hellinger_eculidean.pdf")

plot(ord$points[,1], ord$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat$col, 0.8),
     pch = 16,
     main = "",
     ylim = c(-0.8, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.8,0.4,0.2), labels =  seq(-0.8,0.4,0.2), las = 2)
text(x = -0.5, y=0.53, "a)", xpd = NA, cex = 1.7)

legend(x = 0.4, y = 0.5, 
       legend = unique(merged_dat$habit),
       col = unique(add.alpha(merged_dat$col, 0.8)),
       pch = 19,
       #horiz = T,
       cex = 2,
       xpd = NA,
       bty = "n")
dev.off()

# ##############
#bacteria
#########

rm(list = ls())
library(ecodist)
library(vegan)
library(dendextend)

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#################################################################
# Do ordinations for compartment and host taxon for all samples #
#################################################################


dat <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
dat$sample <- gsub("X","", dat$sample)
metadat <- read.csv("./processedData/16Smetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)

#match order
merged_dat <- merge(metadat, dat, by.x = "sample", by.y = "sample", all.y = T)
dim(merged_dat)
dim(dat)
dim(metadat)
tail(names(merged_dat))
head(names(merged_dat),205)
unique(merged_dat$lifehistory)
names(merged_dat)[length(merged_dat)-4]

datr <-  vegan::decostand(merged_dat[,202:(length(merged_dat)-4)], method = "log")

dat_e <- vegan::vegdist(datr, method =  "altGower")

ord <- cmdscale(dat_e, eig = TRUE, k = 2)

## Assign colors for ordination
library(wesanderson)
colors <- wes_palette("FantasticFox1", n =  4, type = "discrete")

merged_dat$col <- colors[1]
merged_dat$col[merged_dat$compartment == "EP"] <- colors[3]

## See if the ordination correlates with total abundance or the abundance of the duds
## the latter where the unidentifiable sequences
cor.test(ord$points[,2], rowSums(merged_dat[,202:(length(merged_dat)-4)]))
cor.test(ord$points[,2], merged_dat$duds)

###################
# ITS EN vs EP ASVs  #
###################

pdf(width = 8, height = 8, file = "./visuals/ordination_16S_enVsep_asvs_gower_counts_NOISD.pdf")

plot(ord$points[,1], ord$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat$col, 0.8),
     pch = 16,
     main = "",
     xlim = c(-0.05,0.06),
     ylim = c(-0.05,0.05),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.05,0.05,0.025), labels =  seq(-0.05,0.05,0.025), las = 2)
text(x = -0.05, y=0.05, "a)", xpd = NA, cex = 1.7)

#ordihull(ord, merged_dat$col, col = unique(merged_dat$col))
legend(x = 0.03, y = 0.05, 
       legend = c("EN","EP"), 
       col = unique(add.alpha(merged_dat$col, 0.8)), 
       pch = 19, 
       horiz = T,
       cex = 2, 
       xpd = NA,
       bty = "n")

dev.off()

colors <- wes_palette("FantasticFox1", n =  8, type = "continuous")
colors <- colors[c(1,3,4,7)]

k <- 1
for(i in unique(merged_dat$habit)){
  merged_dat$col[merged_dat$habit == i] <- colors[k]
  k <- k + 1
}

pdf(width = 8, height = 8, file = "./visuals/ordination_16S_habit_asvs_gower_counts_NOISD.pdf")

plot(ord$points[,1], ord$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat$col, 0.8),
     pch = 16,
     main = "",
     xlim = c(-0.05,0.06),
     ylim = c(-0.05,0.05),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.05,0.05,0.025), labels =  seq(-0.05,0.05,0.025), las = 2)
text(x = -0.05, y=0.05, "a)", xpd = NA, cex = 1.7)

#ordihull(ord, merged_dat$col, col = unique(merged_dat$col))
legend(x = 0.03, y = -0.01, 
       legend = unique(merged_dat$habit), 
       col = unique(add.alpha(merged_dat$col, 0.8)), 
       pch = 19, 
       #horiz = T,
       cex = 2, 
       xpd = NA,
       bty = "n")

dev.off()

##BRING IN THE ISD

merged_dat[,202:(length(merged_dat)-4)] <- merged_dat[,202:(length(merged_dat)-4)]/merged_dat$ISD

datr <-  vegan::decostand(merged_dat[,202:(length(merged_dat)-4)], method = "hellinger")

dat_e <- vegan::vegdist(datr, method =  "euclidean")

ord <- cmdscale(dat_e, eig = TRUE, k = 2)

cor.test(ord$points[,2], rowSums(merged_dat[,202:(length(merged_dat)-4)]))
cor.test(ord$points[,1], merged_dat$ISD)
cor.test(ord$points[,1], merged_dat$duds)
cor.test(ord$points[,2], merged_dat$ISD + merged_dat$duds)

pdf(width = 8, height = 8, file = "./visuals/ordination_16S_habit_asvs_gower_counts_ISD.pdf")

plot(ord$points[,1], ord$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat$col, 0.8),
     pch = 16,
     main = "",
    # xlim = c(-0.05,0.05),
    # ylim = c(-0.05,0.05),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n"
    )
axis(side = 2, at = c(-0.3,-0.2,-0.1,0,.1), labels =   c(-0.3,-0.2,-0.1,0,.1), las = 2)
text(x = -1.2, y=0.14, "a)", xpd = NA, cex = 1.7)

#ordihull(ord, merged_dat$col, col = unique(merged_dat$col))
legend(x = -1.1, y = 0.13, 
       legend = c("EN","EP"), 
       col = unique(add.alpha(merged_dat$col, 0.8)), 
       pch = 19, 
       horiz = T,
       cex = 2, 
       xpd = NA,
       bty = "n")

dev.off()

######################################
# 
# ## Redo the colors so that they are by host taxon
# #colors <- randomcoloR::randomColor(count = length(unique(merged_dat$taxon_final)))
# colors <- rainbow(n = length(unique(merged_dat$taxon_final)))
# k <- 1
# for(i in unique(merged_dat$taxon_final)){
#   merged_dat$col[merged_dat$taxon_final == i] <- colors[k]
#   k <- k + 1
# }
# 
# ###################
# # ITS by host taxon - ASVs  #
# ###################
# 
# pdf(width = 8, height = 8, file = "./visuals/ordination_16S_byHost_asvshellinger.pdf")
# 
# plot(ord$points[,1], ord$points[,2], 
#      #type = "n", 
#      col = add.alpha(merged_dat$col, 0.5),
#      pch = 16,
#      main = "",
#      ylim = c(-0.7, 0.4),
#      cex = 1.3,
#      frame.plot = F,
#      xlab = "PCoA1", ylab = "PCoA2",
#      yaxt = "n")
# axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
# text(x = -0.1, y=0.6, "b)", xpd = NA, cex = 1.7)
# 
# #Ordihull is not working correctly. Do not feel like debugging.
# #ordihull(ord, merged_dat$taxon_final, col = unique(merged_dat$col))
# dev.off()
# 
# ###################
# # Legend plot #
# ###################
# 
# pdf(width = 7, height = 14, file = "./visuals/ordination16s_LEGENDhellinger.pdf")
# 
# plot(ord$points[,1], ord$points[,2], type = "n",axes = F, ylab = "", xlab = "")
# legend(x = 0, y = 0.2, 
#        legend = unique(merged_dat$taxon_final), 
#        col = unique(add.alpha(merged_dat$col, 0.5)), 
#        pch = 19, 
#        cex = 0.6, 
#        xpd = NA,
#        bty = "n")
# 
# dev.off()
# 
# ###################
# # ITS by life history  #
# ###################
# 
# #colors <- rainbow(n = length(unique(merged_dat$lifehistory)))
# library(wesanderson)
# colors <- wes_palette("FantasticFox1", n =  8, type = "continuous")
# colors <- colors[c(1,3,4,7)]
# 
# k <- 1
# for(i in unique(merged_dat$lifehistory)){
#   merged_dat$col[merged_dat$lifehistory == i] <- colors[k]
#   k <- k + 1
# }
# 
# pdf(width = 8, height = 8, file = "./visuals/ordination_16S_lifehistory_asvshellinger.pdf")
# 
# plot(ord$points[,1], ord$points[,2], 
#      #type = "n", 
#      col = add.alpha(merged_dat$col, 0.8),
#      pch = 16,
#      main = "",
#      ylim = c(-0.7, 0.4),
#      cex = 1.3,
#      frame.plot = F,
#      xlab = "PCoA1", ylab = "PCoA2",
#      yaxt = "n")
# axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
# text(x = -0.1, y=0.6, "b)", xpd = NA, cex = 1.7)
# 
# #Ordihull is not working correctly. Do not feel like debugging.
# #ordihull(ord, merged_dat$taxon_final, col = unique(merged_dat$col))
# dev.off()
# 
# ###################
# # Legend plot #
# ###################
# 
# pdf(width = 8, height = 8, file = "./visuals/ordination16sLIFEHISTORY_LEGENDhellinger.pdf")
# 
# plot(ord$points[,1], ord$points[,2], type = "n",axes = F, ylab = "", xlab = "")
# legend(x = 0, y = 0.2, 
#        legend = unique(merged_dat$lifehistory), 
#        col = unique(add.alpha(merged_dat$col, 0.5)), 
#        pch = 19, 
#        cex = 0.6, 
#        xpd = NA,
#        bty = "n")
# 
# dev.off()
# # ##############
