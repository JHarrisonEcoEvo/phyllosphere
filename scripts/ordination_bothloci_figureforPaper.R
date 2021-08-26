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

datr <-  vegan::decostand(merged_dat[,206:(length(merged_dat)-2)], method = "hellinger")

## Experimentation with other distance metrics and ordination techniques
# dat_bc <- bcdist(datr)
# ord <- pco(x = dat_bc)
#ord$vectors[,1] > 0.001
# plot(ord$vectors[,1], ord$vectors[,2], 
#      #type = "n", 
#      col = metadat$col,
#      xlab = "PCoA1", ylab = "PCoA2",
#      axes = TRUE, main = "PCoA (ecodist)")

#dat_bc <- bcdist(merged_dat[,204:(length(merged_dat)-2)])

dat_e <- dist(datr,method = "euclidean")

# ord <- decorana(dat_e)
#ord <- pco(dat_e)
# plot(ord$vectors[,1], ord$vectors[,2], 
#      #type = "n", 
#      col = merged_dat$col,
#      pch = 16, 
#      xlab = "PCoA1", ylab = "PCoA2",
#      axes = TRUE, main = "PCoA (ecodist)")

ord <- cmdscale(dat_e, eig = TRUE, k = 2)


## Assign colors for ordination
merged_dat$col <- "blue"
merged_dat$col[merged_dat$compartment == "EP"] <-"orange"

## See if the ordination correlates with total abundance or the abundance of the duds
## the latter where the unidentifiable sequences
cor.test(ord$points[,2], rowSums(merged_dat[,206:(length(merged_dat)-2)]))
cor.test(ord$points[,2], merged_dat$duds)

###################
# ITS EN vs EP ASVs  #
###################

pdf(width = 7, height = 7, file = "./visuals/ordination_ITS_enVsep_asvs.pdf")

plot(ord$points[,1], ord$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat$col, 0.5),
     pch = 16,
     main = "",
     ylim = c(-0.7, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
text(x = -0.1, y=0.6, "a)", xpd = NA, cex = 1.7)

#ordihull(ord, merged_dat$col, col = unique(merged_dat$col))
legend(x = 0.1, y = 0.7, 
       legend = c("EN","EP"), 
       col = unique(add.alpha(merged_dat$col, 0.5)), 
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

pdf(width = 7, height = 7, file = "./visuals/ordination_ITS_byHost_asvs.pdf")

plot(ord$points[,1], ord$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat$col, 0.5),
     pch = 16,
     main = "",
     ylim = c(-0.7, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
text(x = -0.1, y=0.6, "b)", xpd = NA, cex = 1.7)

#Ordihull is not working correctly. Do not feel like debugging.
#ordihull(ord, merged_dat$taxon_final, col = unique(merged_dat$col))
dev.off()

###################
# Legend plot #
###################

pdf(width = 7, height = 14, file = "./visuals/ordination_LEGEND.pdf")

plot(ord$points[,1], ord$points[,2], type = "n",axes = F, ylab = "", xlab = "")
legend(x = 0, y = 0.2, 
       legend = unique(merged_dat$taxon_final), 
       col = unique(add.alpha(merged_dat$col, 0.5)), 
       pch = 19, 
       cex = 0.6, 
       xpd = NA,
       bty = "n")

dev.off()

################################################################
# ITS by host taxon - summing ASVs and making a dendrogram  #
################################################################

summed_df <- data.frame(matrix(nrow = 62,ncol = 1))
k <- 1
for(i in 206:(length(merged_dat)-2)){
  summed_df[,k] <- aggregate(merged_dat[,i] ~ merged_dat$taxon_final, FUN = sum)[,2]
  k<- k + 1
}

datr_summed <-  vegan::decostand(summed_df, method = "hellinger")
dat_e_summed <- dist(datr_summed,method = "euclidean")

a <- hclust(d = dat_e_summed,
            method =  "ward.D")

#recall that the tip order is stored in this:
#a$order, which we will use to assign meaningful labels. 

dend1 <- as.dendrogram(a)

dend1 <- dend1 %>% 
  color_labels(col = unique(add.alpha(merged_dat$col, 1))[a$order]) %>%
  set("labels_cex",0.5)

labels(dend1) <- unique(merged_dat$taxon_final)[a$order]

pdf(width = 7, height = 7, file = "./visuals/ordination_ITS_summedByHost_dendrogram.pdf")
par(mar=c(9,3,3,3))

plot(dend1, xpd = NA)

dev.off()


################################################################
# ITS by host taxon - summed ASVs Ordination  #
################################################################

pdf(width = 7, height = 7, file = "./visuals/ordination_ITS_summedByHost_taxon.pdf")


ord <- cmdscale(dat_e_summed, eig = TRUE, k = 2)

plot(ord$points[,1], ord$points[,2],
     #type = "n", 
     col = unique(add.alpha(merged_dat$col, 0.5)),
     pch = 16,
     main = "",
     ylim = c(-0.7, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
text(x = -0.1, y=0.6, "b)", xpd = NA, cex = 1.7)


dev.off()

###########################
# bacteria data wrangling #
###########################

dat16s <- read.csv("./processedData/16sp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv",
                stringsAsFactors = F)

metadat16s <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)

#match order
merged_dat16s <- merge(metadat16s, dat16s, by.x = "sample", by.y = "sample", all.y = T)
dim(merged_dat16s)

datr16s <-  vegan::decostand(merged_dat16s[,204:(length(merged_dat16s)-4)], method = "hellinger")

dat_e16s <- dist(datr16s,method = "euclidean")

ord16s <- cmdscale(dat_e16s, eig = TRUE, k = 2)

merged_dat16s$col <- "blue"
merged_dat16s$col[merged_dat16s$compartment == "EP"] <-"orange"

###########################
# bacteria EP vs EN ordination #
###########################

pdf(width = 7, height = 7, file = "./visuals/ordination_16S_epVsen_asvs.pdf")

plot(ord16s$points[,1], ord16s$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat16s$col, 0.5),
     pch = 16,
     main = "",
     ylim = c(-0.7, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
text(x = -0.6, y=0.5, "e)", xpd = NA, cex = 1.7)

#ordihull(ord, merged_dat$col, col = unique(merged_dat$col))
legend("top",
       legend = c("EN","EP"), 
       col = unique(add.alpha(merged_dat16s$col, 0.5)), 
       pch = 19, 
       horiz = T,
       cex = 2, 
       xpd = NA,
       bty = "n")

dev.off()


###########################
# bacteria EP vs EN ordination zoomed in #
###########################

pdf(width = 7, height = 7, file = "./visuals/ordination_16S_epVsen_asvs_zoomedIn.pdf")

plot(ord16s$points[ord16s$points[,1] > -0.02,1], ord16s$points[ord16s$points[,1] > -0.02,2], 
     #type = "n", 
     col = add.alpha(merged_dat16s$col, 0.5),
     pch = 16,
     main = "",
     ylim = c(-0.7, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
text(x = -0.1, y=0.6, "e)", xpd = NA, cex = 1.7)

#ordihull(ord, merged_dat$col, col = unique(merged_dat$col))
legend(x = -0.4, y = 0.7, 
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

colors <- rainbow(n = length(unique(merged_dat16s$taxon_final)))
k <- 1
for(i in unique(merged_dat16s$taxon_final)){
  merged_dat16s$col[merged_dat16s$taxon_final == i] <- colors[k]
  k <- k + 1
}

###########################
# bacteria by host - asvs #
###########################

pdf(width = 7, height = 7, file = "./visuals/ordination_16S_host_asvs.pdf")

plot(ord16s$points[,1], ord16s$points[,2], 
     #type = "n", 
     col = add.alpha(merged_dat16s$col, 0.5),
     pch = 16,
     main = "",
     ylim = c(-0.7, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
text(x = -0.1, y=0.6, "f)", xpd = NA, cex = 1.7)

dev.off()

#Zoomed in
pdf(width = 7, height = 7, file = "./visuals/ordination_16S_host_asvsZoomedin.pdf")

plot(ord16s$points[ord16s$points[,1] > -0.04,1], ord16s$points[ord16s$points[,1] > -0.04,2], 
     #type = "n", 
     col = add.alpha(merged_dat16s$col, 0.5),
     pch = 16,
     main = "",
     ylim = c(-0.7, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
text(x = -0.1, y=0.6, "f)", xpd = NA, cex = 1.7)

dev.off()

#Lump by host
summed_df16 <- data.frame(matrix(nrow = 63,ncol = 1))
k <- 1
for(i in 204:(length(merged_dat16s)-4)){
  summed_df16[,k] <- aggregate(merged_dat16s[,i] ~ merged_dat16s$taxon_final, FUN = sum)[,2]
  k<- k + 1
}

datr <-  vegan::decostand(summed_df16, method = "hellinger")
dat_e16 <- dist(datr,method = "euclidean")

a <- hclust(d = dat_e16,
            method =  "ward.D")

dend <- as.dendrogram(a)

dend <- dend %>% 
  color_labels(col = unique(add.alpha(merged_dat$col, 1))[a$order]) %>%
  set("labels_cex",0.5)
  
labels(dend) <- unique(merged_dat$taxon_final)[a$order]

###########################
# bacteria by host - asvs - dendrogram #
###########################
pdf(width = 7, height = 7, file = "./visuals/ordination_16S_summedByHost_dendrogram.pdf")

par(mar=c(9,3,3,3))
plot(dend, xpd = NA)
text(x = -0.1, y=0.6, "g)", xpd = NA, cex = 1.7)

dev.off()

###########################
# bacteria by host - summed asvs - ordination #
###########################

pdf(width = 7, height = 7, file = "./visuals/ordination_ITS_summedByHost_taxon.pdf")


ord <- cmdscale(dat_e16, eig = TRUE, k = 2)

plot(ord$points[,1], ord$points[,2],
     #type = "n", 
     col = unique(add.alpha(merged_dat$col, 0.5)),
     pch = 16,
     main = "",
     ylim = c(-0.7, 0.4),
     cex = 1.3,
     frame.plot = F,
     xlab = "PCoA1", ylab = "PCoA2",
     yaxt = "n")
axis(side = 2, at = seq(-0.6,0.4,0.2), labels =  seq(-0.6,0.4,0.2))
text(x = -0.1, y=0.6, "b)", xpd = NA, cex = 1.7)


dev.off()

