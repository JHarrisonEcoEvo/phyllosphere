
#############
# the fungi #
#############

rm(list=ls())
library(tidyverse)
library(plyr)
library("wesanderson")

its <- read.table("processedData/smallmem97_ITS.sintax", fill = T)
head(its)

its_otus <- read.csv("./processedData/ITSp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv",
                stringsAsFactors = F)
length(names(its_otus)) -4

its_meta <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)

merged_dat <- merge(its_meta, its_otus, by.x = "sample", by.y = "sample", all.y = T)

#recode the names of merged dat for each Zotu as the taxonomy hypothesis
for(i in 205:length(merged_dat)){
  if(any(grepl(names(merged_dat)[i], its$V1))){
    names(merged_dat)[i] <- paste(names(merged_dat)[i],
                                  its$V4[its$V1 == names(merged_dat)[i]],
                                  sep = "_")
  }
}

#Extract unique phyla and classes
phyla <- gsub("d:Fungi,p:(\\w+),.*","\\1",its$V4)
phyla[phyla == "d:Fungi"] <- "Unknown"
phyla[phyla == ""] <- "Unknown"
phyla <- gsub("d:Fungi,p:(\\w+)","\\1",phyla)

#Sum abundances by compartment and phylum
phyla_df <- list()
k <- 1
for(i in unique(phyla)){
  if(length(grep(i, names(merged_dat))) > 1){
    df <- aggregate(
    rowSums(merged_dat[, grep(i, names(merged_dat))]) ~
      merged_dat$compartment, FUN = mean)
    df$taxon <- i
    names(df) <- c("compartment","sum","taxon")
    phyla_df[[k]] <- df
  }else if(length(grep(i, names(merged_dat))) == 1){
    df <- aggregate(
    merged_dat[, grep(i, names(merged_dat))] ~
      merged_dat$compartment, FUN = mean)
    df$taxon <- i
    names(df) <- c("compartment","sum","taxon")
    phyla_df[[k]] <- df
  }
  k <- k + 1
}

phyla_df <- phyla_df[-2]

df <- ldply(phyla_df, data.frame)

wide = df %>% 
  spread(compartment, sum)
wide
fornext <- wide$taxon ######### This is to ensure the colors match between plots

colors <- wes_palette("FantasticFox1", n = dim(wide)[1], type = "continuous")

pdf(width = 5,height = 6, file = "./visuals/phyla_its_compartment_barplot_modeledData.pdf")
# Get the stacked barplot
par(mar = c(3,6,1,1), oma = c(2,2,2,2))
barplot(as.matrix(wide[,2:3]), 
        col=colors, 
        border="white", 
        cex.lab = 1.5,
        cex.axis = 1.5,
        las = 2, 
       # ylim = c(0,5),
        space=0.04, 
        font.axis=2, 
        xaxt = "n",
        xlab="")
text( c("EN", "EP"), cex = 1.5, x = c(0.5,1.5), y = -0.5, xpd = NA)
title(ylab = "Abundance (ratio with ISD)", line = 5, cex.lab = 2, xpd = NA)
dev.off()

####################
# Do for host taxa #
####################

#remove the unknowns
merged_dat <- merged_dat[-grep("Unknown",merged_dat$taxon_final),]

#Sum abundances by c
phyla_df <- list()
k <- 1
for(i in unique(phyla)){
  if(length(grep(i, names(merged_dat))) > 1){
    df <- aggregate(
      rowSums(merged_dat[, grep(i, names(merged_dat))]) ~
        merged_dat$taxon_final + merged_dat$compartment, FUN = mean)
    df$taxon <- i
    names(df) <- c("host","compartment", "sum","taxon")
    phyla_df[[k]] <- df
  }else if(length(grep(i, names(merged_dat))) == 1){
    df <- aggregate(
      merged_dat[, grep(i, names(merged_dat))] ~
        merged_dat$taxon_final + merged_dat$compartment, FUN = mean)
    df$taxon <- i
    names(df) <- c("host","compartment", "sum","taxon")
    phyla_df[[k]] <- df
  }
  k <- k + 1
}


phyla_df <- phyla_df[-2]

df <- ldply (phyla_df, data.frame)
df$group <- paste(df$host, df$compartment, sep = "")

library(data.table)
setDT(df)   # coerce to data.table
data_wide <- dcast(df, group ~ taxon, 
                   value.var = c("sum"))

names(data_wide)[-1] == fornext
data_wide_t <- data.frame(t(data_wide[,2:length(data_wide)]))
names(data_wide_t) <- data_wide$group

colors <- wes_palette("FantasticFox1", n = dim(wide)[1], type = "continuous")

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

pdf(width = 12,height = 9, file = "./visuals/phyla_its_taxon_final_barplot_modeledData.pdf")
# Get the stacked barplot
par(mar = c(3,6,1,1), oma = c(16,3,1,1))
a <- barplot(as.matrix(data_wide_t), 
        col = colors,
        border="white", 
        cex.lab = 1.5,
        cex.axis = 1.5,
        las = 2, 
        ylim = c(0,100),
        space=0.04, 
        #font.axis=2, 
        xaxt = "n",
        ylab = "",
        xlab="")

k <- 1
a_means <- NA
for(i in seq(1,length(names(data_wide_t)), by = 2)){
  a_means[k] <- mean(a[i:(i+1)])
  k <- k + 1
}

rect(a_means, 
     0, 
     1 + a_means, 
     colSums(data_wide_t)[seq(2,length(data_wide_t), 2)], 
     density = 50, angle = 45,
     col = alpha("black", 0.4),
     border = NA)

text(gsub("E[NP]", "",names(data_wide_t))[
  seq(1,length(names(data_wide_t)), by = 2)],
     pos = 2,
     srt = 90,
     cex = 0.5, 
     font = 3,
     x = 1.3 + a_means,
     y = -0.1,
     xpd = NA)
title(ylab = "Abundance (ratio with ISD)", line = 4, cex.lab = 2)
#create positions for tick marks, one more than number of bars
dev.off()

##########
# legend #
##########
pdf(width = 3,height = 6, file = "./visuals/phyla_its_barplot_legend_modeledData.pdf")

plot(x = seq(1,6,1), y = seq(1,6,1),
     type = "n",
     frame.plot = F,
     axes = F,
     xlab = "",
     ylab = "")
legend(legend = wide$taxon,
       x = "left",
       col = colors,
       xpd = NA,
       bty = "n",
       pch = 15)
dev.off()

####################################################################

#############
# the bacteria #
#############

rm(list=ls())

bacts <- read.table("processedData/smallmem97_16S.sintax", fill = T)

bact_otus <- read.csv("./processedData/16sp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv",
                     stringsAsFactors = F)

bact_meta <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                     stringsAsFactors = F)

merged_dat <- merge(bact_meta, bact_otus, by.x = "sample", by.y = "sample", all.y = T)

#recode the names of merged dat for each Zotu as the taxonomy hypothesis
for(i in 204:length(merged_dat)){
  if(any(grepl(names(merged_dat)[i], bacts$V1))){
    names(merged_dat)[i] <- paste(names(merged_dat)[i],
                                  bacts$V4[bacts$V1 == names(merged_dat)[i]],
                                  sep = "_")
  }
}

#Extract unique phyla and classes
phyla <- gsub("k:Bacteria,p:(\\w+),.*","\\1",bacts$V4)
phyla[phyla == "k:Bacteria"] <- "Unknown"
phyla[phyla == ""] <- "Unknown"
phyla <- gsub("k:Bacteria,p:(\\w+)","\\1",phyla)
phyla <- gsub("k:Archaea,p:(\\w+),.*","Archaea-\\1",phyla)
phyla <- gsub("k:Archaea,p:(\\w+)","Archaea-\\1",phyla)
phyla <- gsub("k:Bacteria,p:\\[(\\w+)\\],.*","\\1",phyla)
phyla <- gsub("k:Archaea,p:\\[(\\w+)\\],.*","\\1",phyla)

unique(phyla)

####################
# Do for host taxa #
####################

#remove the unknowns
merged_dat <- merged_dat[-grep("Unknown",merged_dat$taxon_final),]

#Sum abundances by c
phyla_df <- list()
k <- 1
for(i in unique(phyla)){
  if(length(grep(i, names(merged_dat))) > 1){
    df <- aggregate(
      rowSums(merged_dat[, grep(i, names(merged_dat))]) ~
        merged_dat$taxon_final + merged_dat$compartment, FUN = mean)
    df$taxon <- i
    names(df) <- c("host","compartment", "sum","taxon")
    phyla_df[[k]] <- df
  }else if(length(grep(i, names(merged_dat))) == 1){
    df <- aggregate(
      merged_dat[, grep(i, names(merged_dat))] ~
        merged_dat$taxon_final + merged_dat$compartment, FUN = mean)
    df$taxon <- i
    names(df) <- c("host","compartment", "sum","taxon")
    phyla_df[[k]] <- df
  }
  k <- k + 1
}


phyla_df <- phyla_df[-2]

df <- ldply (phyla_df, data.frame)
df$group <- paste(df$host, df$compartment, sep = "")

library(data.table)
setDT(df)   # coerce to data.table
data_wide <- dcast(df, group ~ taxon, 
                   value.var = c("sum"))

data_wide_t <- data.frame(t(data_wide[,2:length(data_wide)]))
names(data_wide_t) <- data_wide$group

library("wesanderson")
colors <- wes_palette("FantasticFox1", n = dim(data_wide_t)[1], type = "continuous")
#rearrange colors so that the abundant stuff have differentiable colors
wide <- data_wide_t[order(rowSums(data_wide_t[2:length(data_wide_t)])),]
fornext <- row.names(wide)

hidef <- seq(1,47,5)
colors_new <- rep(colors[hidef[10]], 47)
colors_new[39:47] <- colors[hidef[1:9]]

pdf(width = 12,height = 9, file = "./visuals/phyla_16s_taxon_final_barplot_modeledData.pdf")
# Get the stacked barplot
par(mar = c(3,6,1,1), oma = c(16,3,1,1))

#Sanity check, plug this in to see stuff change
#wide[47,10] <- 20000
a <- barplot(as.matrix(wide), 
             col=colors_new, 
             border="white", 
             cex.lab = 1.5,
             cex.axis = 1.5,
             las = 2, 
             ylim = c(0,300),
             space=0.04, 
             #font.axis=2, 
             xaxt = "n",
             ylab = "",
             xlab="")


k <- 1
a_means <- NA
for(i in seq(1,length(names(wide)), by = 2)){
  a_means[k] <- mean(a[i:(i+1)])
  k <- k + 1
}

rect(a_means, 
     5, 
     1 + a_means, 
     colSums(wide)[seq(2,length(wide), 2)], 
     density = 50, angle = 45,
     col = alpha("black", 0.4),
     border = NA)

text(gsub("E[NP]", "",names(data_wide_t))[
  seq(1,length(names(data_wide_t)), by = 2)],
  pos = 2,
  srt = 90,
  cex = 0.5, 
  font = 3,
  x = 1.3 + a_means,
  y = -0.1,
  xpd = NA)

title(ylab = "Abundance (ratio with ISD)", line = 5, cex.lab = 2, xpd = NA)
#create positions for tick marks, one more than number of bars
dev.off()

#####################
# Do by compartment #
#####################

#reloading this bc I removed the unknowns above
merged_dat <- merge(bact_meta, bact_otus, by.x = "sample", by.y = "sample", all.y = T)

#recode the names of merged dat for each Zotu as the taxonomy hypothesis
for(i in 204:length(merged_dat)){
  if(any(grepl(names(merged_dat)[i], bacts$V1))){
    names(merged_dat)[i] <- paste(names(merged_dat)[i],
                                  bacts$V4[bacts$V1 == names(merged_dat)[i]],
                                  sep = "_")
  }
}

#Sum abundances by compartment and phylum
phyla_df <- list()
k <- 1
for(i in unique(phyla)){
  if(length(grep(i, names(merged_dat))) > 1){
    df <- aggregate(
      rowSums(merged_dat[, grep(i, names(merged_dat))]) ~
        merged_dat$compartment, FUN = mean)
    df$taxon <- i
    names(df) <- c("compartment","sum","taxon")
    phyla_df[[k]] <- df
  }else if(length(grep(i, names(merged_dat))) == 1){
    df <- aggregate(
      merged_dat[, grep(i, names(merged_dat))] ~
        merged_dat$compartment, FUN = mean)
    df$taxon <- i
    names(df) <- c("compartment","sum","taxon")
    phyla_df[[k]] <- df
  }
  k <- k + 1
}

phyla_df <- phyla_df[-2]

df <- ldply (phyla_df, data.frame)

library(tidyverse)
wide = df %>% 
  spread(compartment, sum)
wide
wide <- wide[match( fornext,wide$taxon),]

pdf(width = 5,height = 6, file = "./visuals/phyla_16s_compartment_barplot_modeledData.pdf")
# Get the stacked barplot
par(mar = c(3,6,1,1), oma = c(2,2,2,2))
barplot(as.matrix(wide[,2:3]), 
        col=colors_new, 
        border="white", 
        cex.lab = 1.5,
        cex.axis = 1.5,
        las = 2, 
        ylim = c(0,150),
        space=0.04, 
        font.axis=2, 
        xaxt = "n",
        xlab="")
text( c("EN", "EP"), cex = 1.5, x = c(0.5,1.5), y = -4, xpd = NA)
title(ylab = "Abundance (ratio with ISD)", line = 5, cex.lab = 2, xpd = NA)
dev.off()


##########
# legend #
##########
pdf(width = 3,height = 14, file = "./visuals/phyla_16s_barplot_legend_modeledData.pdf")
par(oma = c(10,2,10,1))
plot(x = seq(1,10,1), y = seq(1,10,1),
     type = "n",
     frame.plot = F,
     axes = F,
     xlab = "",
     ylab = "")
legend(legend = rev(wide$taxon),
       x = "left",
       col = rev(colors_new),
       xpd = NA,
       bty = "n",
       pch = 15)
dev.off()

#################################################
# Do with unnormalized data 
#################################################
#################################################
#################################################
#################################################
#################################################
# #############
# # the fungi #
# #############
# 
# rm(list=ls())
# 
# its <- read.table("processedData/smallmem97_ITS.sintax", fill = T)
# head(its)
# 
# its_otus <- read.csv("./processedData/ITSp_estimates_wrangled_for_post_modeling_analysis.csv",
#                      stringsAsFactors = F)
# 
# its_meta <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
#                      stringsAsFactors = F)
# 
# merged_dat <- merge(its_meta, its_otus, by.x = "sample", by.y = "sample", all.y = T)
# 
# #recode the names of merged dat for each Zotu as the taxonomy hypothesis
# for(i in 205:length(merged_dat)){
#   if(any(grepl(names(merged_dat)[i], its$V1))){
#     names(merged_dat)[i] <- paste(names(merged_dat)[i],
#                                   its$V4[its$V1 == names(merged_dat)[i]],
#                                   sep = "_")
#   }
# }
# 
# #Extract unique phyla and classes
# phyla <- gsub("d:Fungi,p:(\\w+),.*","\\1",its$V4)
# phyla[phyla == "d:Fungi"] <- "Unknown"
# phyla[phyla == ""] <- "Unknown"
# phyla <- gsub("d:Fungi,p:(\\w+)","\\1",phyla)
# 
# #Sum abundances by compartment and phylum
# phyla_df <- list()
# k <- 1
# for(i in unique(phyla)){
#   if(length(grep(i, names(merged_dat))) > 1){
#     df <- aggregate(
#       rowSums(merged_dat[, grep(i, names(merged_dat))]) ~
#         merged_dat$compartment, FUN = mean)
#     df$taxon <- i
#     names(df) <- c("compartment","sum","taxon")
#     phyla_df[[k]] <- df
#   }else if(length(grep(i, names(merged_dat))) == 1){
#     df <- aggregate(
#       merged_dat[, grep(i, names(merged_dat))] ~
#         merged_dat$compartment, FUN = mean)
#     df$taxon <- i
#     names(df) <- c("compartment","sum","taxon")
#     phyla_df[[k]] <- df
#   }
#   k <- k + 1
# }
# 
# phyla_df <- phyla_df[-2]
# 
# df <- ldply (phyla_df, data.frame)
# 
# library(tidyverse)
# wide = df %>% 
#   spread(compartment, sum)
# wide
# fornext <- wide$taxon
# 
# library("wesanderson")
# colors <- wes_palette("FantasticFox1", n = dim(wide)[1], type = "continuous")
# 
# pdf(width = 5,height = 6, file = "./visuals/phyla_its_compartment_barplot_unnormalized.pdf")
# # Get the stacked barplot
# par(mar = c(3,6,1,1), oma = c(2,2,2,2))
# barplot(as.matrix(wide[,2:3]), 
#         col=colors, 
#         border="white", 
#         cex.lab = 1.5,
#         cex.axis = 1.5,
#         las = 2, 
#        # ylim = c(0,10000),
#         space=0.04, 
#         font.axis=2, 
#         xaxt = "n",
#         xlab="")
# text( c("EN", "EP"), cex = 1.5, x = c(0.5,1.5), y = -0.8, xpd = NA)
# title(ylab = "Abundance (ratio with ISD)", line = 5, cex.lab = 2, xpd = NA)
# dev.off()
# 
# ####################
# # Do for host taxa #
# ####################
# 
# #remove the unknowns
# merged_dat <- merged_dat[-grep("Unknown",merged_dat$taxon_final),]
# 
# #Sum abundances by compartment and phylum
# phyla_df <- list()
# k <- 1
# for(i in unique(phyla)){
#   if(length(grep(i, names(merged_dat))) > 1){
#     df <- aggregate(
#       rowSums(merged_dat[, grep(i, names(merged_dat))]) ~
#         merged_dat$taxon_final, FUN = mean)
#     df$taxon <- i
#     names(df) <- c("host","sum","taxon")
#     phyla_df[[k]] <- df
#   }else if(length(grep(i, names(merged_dat))) == 1){
#     df <- aggregate(
#       merged_dat[, grep(i, names(merged_dat))] ~
#         merged_dat$taxon_final, FUN = mean)
#     df$taxon <- i
#     names(df) <- c("host","sum","taxon")
#     phyla_df[[k]] <- df
#   }
#   k <- k + 1
# }
# 
# 
# phyla_df <- phyla_df[-2]
# 
# df <- ldply (phyla_df, data.frame)
# 
# library(tidyverse)
# wide = df %>% 
#   spread(host, sum)
# wide
# wide <- wide[match( fornext,wide$taxon),]
# 
# colors <- wes_palette("FantasticFox1", n = dim(wide)[1], type = "continuous")
# 
# pdf(width = 12,height = 9, file = "./visuals/phyla_its_taxon_final_barplot_unnormalized.pdf")
# # Get the stacked barplot
# par(mar = c(3,6,1,1), oma = c(16,3,1,1))
# a <- barplot(as.matrix(wide[,2:length(wide)]), 
#              col=colors, 
#              border="white", 
#              cex.lab = 1.5,
#              cex.axis = 1.5,
#              las = 2, 
#             # ylim = c(0,1000),
#              space=0.04, 
#              #font.axis=2, 
#              xaxt = "n",
#              ylab = "",
#              xlab="")
# text(names(wide)[-1],
#      pos = 2,
#      srt = 90,
#      cex = 1, 
#      font = 3,
#      x = 0.8 + a,
#      y = -0.1,
#      xpd = NA)
# title(ylab = "Abundance (ratio with ISD)", line = 4, cex.lab = 2)
# #create positions for tick marks, one more than number of bars
# dev.off()
# 
# ##########
# # legend #
# ##########
# pdf(width = 3,height = 9, file = "./visuals/phyla_its_barplot_legend_unnormalized.pdf")
# 
# plot(x = seq(1,10,1), y = seq(1,10,1),
#      type = "n",
#      frame.plot = F,
#      axes = F,
#      xlab = "",
#      ylab = "")
# legend(legend = wide$taxon,
#        x = "left",
#        col = colors,
#        xpd = NA,
#        bty = "n",
#        pch = 15)
# dev.off()
# 
# ####################################################################
# 
# #############
# # the bacteria #
# #############
# 
# rm(list=ls())
# 
# bacts <- read.table("processedData/smallmem97_16S.sintax", fill = T)
# 
# bact_otus <- read.csv("./processedData/16sp_estimates_wrangled_for_post_modeling_analysis.csv",
#                       stringsAsFactors = F)
# 
# bact_meta <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
#                       stringsAsFactors = F)
# 
# merged_dat <- merge(bact_meta, bact_otus, by.x = "sample", by.y = "sample", all.y = T)
# 
# #recode the names of merged dat for each Zotu as the taxonomy hypothesis
# for(i in 204:length(merged_dat)){
#   if(any(grepl(names(merged_dat)[i], bacts$V1))){
#     names(merged_dat)[i] <- paste(names(merged_dat)[i],
#                                   bacts$V4[bacts$V1 == names(merged_dat)[i]],
#                                   sep = "_")
#   }
# }
# 
# #Extract unique phyla and classes
# phyla <- gsub("k:Bacteria,p:(\\w+),.*","\\1",bacts$V4)
# phyla[phyla == "k:Bacteria"] <- "Unknown"
# phyla[phyla == ""] <- "Unknown"
# phyla <- gsub("k:Bacteria,p:(\\w+)","\\1",phyla)
# phyla <- gsub("k:Archaea,p:(\\w+),.*","Archaea-\\1",phyla)
# phyla <- gsub("k:Archaea,p:(\\w+)","Archaea-\\1",phyla)
# phyla <- gsub("k:Bacteria,p:\\[(\\w+)\\],.*","\\1",phyla)
# phyla <- gsub("k:Archaea,p:\\[(\\w+)\\],.*","\\1",phyla)
# 
# unique(phyla)
# 
# ####################
# # Do for host taxa #
# ####################
# 
# #remove the unknowns
# merged_dat <- merged_dat[-grep("Unknown",merged_dat$taxon_final),]
# 
# #Sum abundances by compartment and phylum
# phyla_df <- list()
# k <- 1
# for(i in unique(phyla)){
#   if(length(grep(i, names(merged_dat))) > 1){
#     df <- aggregate(
#       rowSums(merged_dat[, grep(i, names(merged_dat))]) ~
#         merged_dat$taxon_final, FUN = mean)
#     df$taxon <- i
#     names(df) <- c("host","sum","taxon")
#     phyla_df[[k]] <- df
#   }else if(length(grep(i, names(merged_dat))) == 1){
#     df <- aggregate(
#       merged_dat[, grep(i, names(merged_dat))] ~
#         merged_dat$taxon_final, FUN = mean)
#     df$taxon <- i
#     names(df) <- c("host","sum","taxon")
#     phyla_df[[k]] <- df
#   }
#   k <- k + 1
# }
# 
# phyla_df <- phyla_df[-2]
# 
# df <- ldply (phyla_df, data.frame)
# 
# wide = df %>% 
#   spread(host, sum)
# wide
# 
# library("wesanderson")
# colors <- wes_palette("FantasticFox1", n = dim(wide)[1], type = "continuous")
# #rearrange colors so that the abundant stuff have differientiable colors
# wide <- wide[order(rowSums(wide[2:length(wide)])),]
# fornext <- wide$taxon
# 
# hidef <- seq(1,47,5)
# colors_new <- rep(colors[hidef[10]], 47)
# colors_new[39:47] <- colors[hidef[1:9]]
# 
# pdf(width = 12,height = 9, file = "./visuals/phyla_16s_taxon_final_barplot_unnorm.pdf")
# # Get the stacked barplot
# par(mar = c(3,6,1,1), oma = c(16,3,1,1))
# 
# #Sanity check
# #wide[47,10] <- 20000
# a <- barplot(as.matrix(wide[,2:length(wide)]), 
#              col=colors_new, 
#              border="white", 
#              cex.lab = 1.5,
#              cex.axis = 1.5,
#              las = 2, 
#              #ylim = c(0,1000),
#              space=0.04, 
#              #font.axis=2, 
#              xaxt = "n",
#              ylab = "",
#              xlab="")
# text(names(wide)[-1],
#      pos = 2,
#      srt = 90,
#      cex = 1, 
#      font = 3,
#      x = 0.8 + a,
#      y = -0.1,
#      xpd = NA)
# title(ylab = "Abundance (ratio with ISD)", line = 5, cex.lab = 2, xpd = NA)
# #create positions for tick marks, one more than number of bars
# dev.off()
# 
# #####################
# # Do by compartment #
# #####################
# 
# #reloading this bc I removed the unknowns above
# merged_dat <- merge(bact_meta, bact_otus, by.x = "sample", by.y = "sample", all.y = T)
# 
# #recode the names of merged dat for each Zotu as the taxonomy hypothesis
# for(i in 204:length(merged_dat)){
#   if(any(grepl(names(merged_dat)[i], bacts$V1))){
#     names(merged_dat)[i] <- paste(names(merged_dat)[i],
#                                   bacts$V4[bacts$V1 == names(merged_dat)[i]],
#                                   sep = "_")
#   }
# }
# 
# #Sum abundances by compartment and phylum
# phyla_df <- list()
# k <- 1
# for(i in unique(phyla)){
#   if(length(grep(i, names(merged_dat))) > 1){
#     df <- aggregate(
#       rowSums(merged_dat[, grep(i, names(merged_dat))]) ~
#         merged_dat$compartment, FUN = mean)
#     df$taxon <- i
#     names(df) <- c("compartment","sum","taxon")
#     phyla_df[[k]] <- df
#   }else if(length(grep(i, names(merged_dat))) == 1){
#     df <- aggregate(
#       merged_dat[, grep(i, names(merged_dat))] ~
#         merged_dat$compartment, FUN = mean)
#     df$taxon <- i
#     names(df) <- c("compartment","sum","taxon")
#     phyla_df[[k]] <- df
#   }
#   k <- k + 1
# }
# 
# phyla_df <- phyla_df[-2]
# 
# df <- ldply (phyla_df, data.frame)
# 
# library(tidyverse)
# wide = df %>% 
#   spread(compartment, sum)
# wide
# wide <- wide[match( fornext,wide$taxon),]
# 
# pdf(width = 5,height = 6, file = "./visuals/phyla_16s_compartment_barplot_unnorm.pdf")
# # Get the stacked barplot
# par(mar = c(3,6,1,1), oma = c(2,2,2,2))
# barplot(as.matrix(wide[,2:3]), 
#         col=colors_new, 
#         border="white", 
#         cex.lab = 1.5,
#         cex.axis = 1.5,
#         las = 2, 
#         #ylim = c(0,10000),
#         space=0.04, 
#         font.axis=2, 
#         xaxt = "n",
#         xlab="")
# text( c("EN", "EP"), cex = 1.5, x = c(0.5,1.5), y = -0.1, xpd = NA)
# title(ylab = "Abundance (ratio with ISD)", line = 5, cex.lab = 2, xpd = NA)
# dev.off()
# 
# 
# ##########
# # legend #
# ##########
# pdf(width = 3,height = 14, file = "./visuals/phyla_16s_barplot_legend_unnorm.pdf")
# par(oma = c(10,2,10,1))
# plot(x = seq(1,10,1), y = seq(1,10,1),
#      type = "n",
#      frame.plot = F,
#      axes = F,
#      xlab = "",
#      ylab = "")
# legend(legend = rev(wide$taxon),
#        x = "left",
#        col = rev(colors_new),
#        xpd = NA,
#        bty = "n",
#        pch = 15)
# dev.off()
