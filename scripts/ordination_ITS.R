rm(list = ls())
library(ecodist)
library(vegan)
dat <- read.csv("./processedData/ITSp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv",
                stringsAsFactors = F)

metadat <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)

#match order
merged_dat <- merge(metadat, dat, by.x = "sample", by.y = "sample", all.y = T)
dim(merged_dat)


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

datr <-  vegan::decostand(merged_dat[,206:(length(merged_dat)-2)], method = "hellinger")

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


#Do plot

layout(matrix(c(1, 2, 3,
                4, 4, 3,
                5,6,3,
                7,7,3), nrow=2, byrow=TRUE))

merged_dat$col <- "blue"
merged_dat$col[merged_dat$compartment == "EP"] <-"orange"

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

merged_dat$col <- "blue"
#colors <- randomcoloR::randomColor(count = length(unique(merged_dat$taxon_final)))
colors <- rainbow(n = length(unique(merged_dat$taxon_final)))
k <- 1
for(i in unique(merged_dat$taxon_final)){
  merged_dat$col[merged_dat$taxon_final == i] <- colors[k]
  k <- k + 1
}

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

#Ordihull is not working
#ordihull(ord, merged_dat$taxon_final, col = unique(merged_dat$col))

plot(ord$points[,1], ord$points[,2], type = "n",axes = F, ylab = "", xlab = "")
legend(x = 0, y = 0.2, 
       legend = unique(merged_dat$taxon_final), 
       col = unique(add.alpha(merged_dat$col, 0.5)), 
       pch = 19, 
       cex = 0.6, 
       xpd = NA,
       bty = "n")

#Lump by host and do a dendrogram
summed_df <- data.frame(matrix(nrow = 62,ncol = 1))
k <- 1
for(i in 206:(length(merged_dat)-2)){
  summed_df[,k] <- aggregate(merged_dat[,i] ~ merged_dat$taxon_final, FUN = sum)[,2]
  k<- k + 1
}

datr <-  vegan::decostand(summed_df, method = "hellinger")
dat_e <- dist(datr,method = "euclidean")


library(dendextend)
#par(mar=c(10,2,2,2), oma =c(2,0,0,0), mfrow=c(1,1))
dend <- as.dendrogram(hclust(d = dat_e,
                             method =  "ward.D"))

dend <- dend %>%
  color_labels(col = unique(add.alpha(merged_dat$col, 1))) 

labels(dend) <- unique(merged_dat$taxon_final)
plot(dend, xpd = NA)

#####################################
# Perform ordinations by host taxon #
#####################################
# 
# merged_dat$col <- "blue"
# merged_dat$col[merged_dat$compartment == "EP"] <-"orange"
# 
# dat_bc <- bcdist(merged_dat[merged_dat$taxon_final == "Artemisia tridentata" &
#                               merged_dat$siteName == "pilot peak",
#                             204:(length(merged_dat)-2)])
# 
# ord <- pco(x = dat_bc)
# plot(ord$vectors[,1], ord$vectors[,2], 
#      #type = "n", 
#      col = merged_dat$col[merged_dat$taxon_final == "Artemisia tridentata"],
#      pch = 16, 
#      xlab = "PCoA1", ylab = "PCoA2",
#      axes = TRUE, main = "PCoA (ecodist)")


############################
#sanity check doing with its90
############################
############################
# 
# dat <- read.table("./processedData/otuTables/its90_otuTable",
#                   stringsAsFactors = F)
# dat[1:3,1:5]
# 
# #Need to get metadata into the same order as the modeled data
# sample <- gsub("_16S_[1,2]_\\w+$", "",dat[1,])
# sample <- gsub("_16S$", "",sample)
# sample <- gsub("dupe_", "",sample)
# sample <- gsub("_rewash$", "",sample)
# sample <- gsub("(E[NP]).*$", "\\1",sample)
# sample <- gsub("(\\d+)(E[NP])$", "\\1_\\2",sample)
# 
# #Need to get metadata into the same order as the modeled data
# sample <- gsub("_ITS_[1,2]_\\w+$", "",dat[1,])
# sample <- gsub("ITS_[ATCG]*_[ATCG]*_", "",sample)
# sample <- toupper(sample)
# sample <- gsub("(\\d+)(E)", "\\1_\\2",sample)
# 
# metadat <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
#                     stringsAsFactors = F)
# 
# table(sample %in% metadat$sample) #some are missing, but good enough for qc
# zotus <- dat[,1]
# 
# dat <- data.frame(t(dat[2:length(dat[,1]),2:length(dat)]))
# dat[,1:3]
# names(dat) <- zotus
# dat$sample <- sample[2:length(sample)]
# 
# metadat$col <- "blue"
# metadat$col[metadat$compartment == "EP"] <-"orange"
# 
# #match order
# merged_dat <- merge(metadat, dat, by.x = "sample", by.y = "sample")
# dim(merged_dat)
# dim(dat)
# dim(metadat)
# 
# #Note there are pcr duplicates. Just ignoring. 
# 
# dat_bc <- ecodist::bcdist(merged_dat[merged_dat$taxon_final == "Artemisia tridentata" &
#                               merged_dat$siteName == "pilot peak",
#                             204:(length(merged_dat)-2)])
# 
# ord <- ecodist::pco(x = dat_bc)
# plot(ord$vectors[,1], ord$vectors[,2], 
#      #type = "n", 
#      col = merged_dat$col[merged_dat$taxon_final == "Artemisia tridentata"],
#      pch = 16, 
#      xlab = "PCoA1", ylab = "PCoA2",
#      axes = TRUE, main = "PCoA (ecodist)")
# 
# #we get a similar pattern
# 
# 
# 
# rm(list = ls())
# library(ecodist)
# dat <- read.csv("./processedData/otuTables/its97_cluster_otus_algorithm_otutableCLEAN",
#                 stringsAsFactors = F, header = T)
# 
# datr <-  vegan::decostand(dat[1:(length(dat$OTUID)-2),2:length(dat)], method = "hellinger")
# 
# datr[1:3,1:5]
# dim(datr)
# 
# dat_bc <- bcdist(datr)
# 
# ord <- pco(x = dat_bc)
# 
# pdf(width = 8, height = 8, file = "its_hellinger_cleanedotutable_ord.pdf")
# plot(ord$vectors[,1], ord$vectors[,2], 
#      xlab = "PCoA1", ylab = "PCoA2",
#      axes = TRUE, main = "PCoA (ecodist)")
# dev.off()
