To do. 
#confirm that the tenericutes are really abundant. What the fuck are they...maybe a grep problems(
#split the axis on the bacteria plot so that primula parryi doesn't outshine all else
#############
# the fungi #
#############

rm(list=ls())
library(tidyverse)
library(plyr)
library("wesanderson")

#SINTAX
its <- read.table("processedData/smallmem97_ITS.sintax", fill = T)
head(its)

its_otus <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
length(names(its_otus)) -4

its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] -1
its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] / its_otus$ISD

its_meta <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)

its_otus$sample <- gsub("X", "", its_otus$sample)

merged_dat <- merge(its_meta, its_otus, by.x = "sample", by.y = "sample", all.y = T)
dim(merged_dat)
dim(its)
dim(its_otus)
names(merged_dat)[200:205]

#recode the names of merged dat for each Zotu as the taxonomy hypothesis
for(i in 204:length(merged_dat)){
  if(any(grepl(names(merged_dat)[i], as.character(its$V1)))){
    names(merged_dat)[i] <- paste(names(merged_dat)[i],
                                  its$V4[its$V1 == names(merged_dat)[i]],
                                  sep = "_")
  }
}

#Extract unique phyla, family, class, genus. 
phyla <- gsub("d:Fungi,p:(\\w+),.*","\\1",its$V4) #Phyla is NOT lowest taxonomy hypothesized
phyla[phyla == "d:Fungi"] <- "Unknown" #recoding unidentified stuff
phyla[phyla == ""] <- "Unknown"
phyla <- gsub("d:Fungi,p:(\\w+)","\\1",phyla) #Phyla is lowest taxonomy hypothesized

class <- gsub("d:Fungi,.*c:(\\w+),.*","\\1",its$V4) #class is NOT lowest taxonomy hypothesized
class[class == "d:Fungi"] <- "Unknown"
class[class == ""] <- "Unknown"
class <- gsub("d:Fungi,.*c:(\\w+)","\\1",class) #class is lowest taxonomy hypothesized
class[grep("p:", class )] <- "Unknown" #these didn't have a class

order <- gsub("d:Fungi,.*o:(\\w+),.*","\\1",its$V4) 
order[order == "d:Fungi"] <- "Unknown"
order[order == ""] <- "Unknown"
order <- gsub("d:Fungi,.*o:(\\w+)","\\1", order) 
order[grep("p:", order )] <- "Unknown" 

family <- gsub("d:Fungi,.*f:(\\w+),.*","\\1",its$V4) 
family[family == "d:Fungi"] <- "Unknown"
family[family == ""] <- "Unknown"
family <- gsub("d:Fungi,.*f:(\\w+)","\\1", family) 
family[grep("p:", family )] <- "Unknown" 

genus <- gsub("d:Fungi,.*g:(\\w+),.*","\\1",its$V4) 
genus[genus == "d:Fungi"] <- "Unknown"
genus[genus == ""] <- "Unknown"
genus <- gsub("d:Fungi,.*g:(\\w+)","\\1", genus) 
genus[grep("p:", genus )] <- "Unknown" 

#Sum abundances by compartment and phylum
taxa_counter <- function(taxa, data, pattern1, pattern2, pattern3) {
  #Pattern 1 and 2 are strings for regex that defines stuff that is unknown.
  outdf <- list()
  k <- 1
  names(data)[grep(pattern1, names(data))] <- "Unknown"
  names(data)[grep(pattern2, names(data))] <- "Unknown"
  names(data)[grep(pattern3, names(data))] <- "Unknown"
  
  for (i in unique(taxa)) {
    if (length(grep(i, names(data))) > 1) {
      df <- aggregate(rowSums(data[, grep(i, names(data))]) ~
                        data$compartment, FUN = median)
      df$taxon <- i
      names(df) <- c("compartment", "sum", "taxon")
      outdf[[k]] <- df
    } else if (length(grep(i, names(data))) == 1) {
      df <- aggregate(data[, grep(i, names(data))] ~
                        data$compartment, FUN = median)
      df$taxon <- i
      names(df) <- c("compartment", "sum", "taxon")
      outdf[[k]] <- df
    }
    k <- k + 1
  }
  return(outdf)
}

phyla_counts <- taxa_counter(taxa = phyla, 
                             data = merged_dat, 
                             pattern1 = ".*d:Fungi$",
                             pattern2 = "Zotu\\d+_$",
                             pattern3 = "Zotu\\d+_$")

class_counts <- taxa_counter(taxa = class, 
                             data = merged_dat, 
                             pattern1 = ".*p:\\w+$",
                             pattern2 = "Zotu\\d+_$",
                             pattern3 = ".*d:Fungi$")
order_counts <- taxa_counter(taxa = order, 
                             data = merged_dat, 
                             pattern1 = ".*p:\\w+$",
                             pattern2 = "Zotu\\d+_$",
                             pattern3 = ".*d:Fungi$")
family_counts <- taxa_counter(taxa = family, 
                             data = merged_dat, 
                             pattern1 = ".*p:\\w+$",
                             pattern2 = "Zotu\\d+_$",
                             pattern3 = ".*d:Fungi$")
genus_counts <- taxa_counter(taxa = genus, 
                             data = merged_dat, 
                             pattern1 = ".*p:\\w+$",
                             pattern2 = "Zotu\\d+_$",
                             pattern3 = ".*d:Fungi$")


#given low counts...I think phyla level plots are the best, with summary stats provided for other notable taxonomic groups

##################

df <- ldply(phyla_counts, data.frame)

wide = df %>% 
  spread(compartment, sum)
wide
fornext <- wide$taxon ######### This is to ensure the colors match between plots

library(wesanderson)
colors <- wes_palette("FantasticFox1", n = dim(wide)[1], type = "continuous")

pdf(width = 5,height = 6, file = "./visuals/phyla_its_compartment_barplot_count.pdf")
# Get the stacked barplot
par(mar = c(3,6,1,1), oma = c(2,2,2,2))
barplot(as.matrix(wide[,2:3]), 
        col=colors, 
        border="white", 
        cex.lab = 1.5,
        cex.axis = 1.5,
        las = 2, 
        ylim = c(0,.1),
        space=0.04, 
        font.axis=2, 
        xaxt = "n",
        xlab="")
text( c("EN", "EP"), cex = 1.5, x = c(0.5,1.5), y = -0.01, xpd = NA)
title(ylab = "Abundance (ratio with ISD)", line = 5, cex.lab = 2, xpd = NA)
dev.off()

#Sanity check to make sure colors don't get out of order
# pdf(width = 3,height = 6, file = "./visuals/phyla_its_barplot_legend_modeledDataTst.pdf")
# 
# plot(x = seq(1,6,1), y = seq(1,6,1),
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
####################
# Do for host taxa #
####################

#remove the unknowns
merged_dat <- merged_dat[-grep("Unknown",merged_dat$taxon_final),]

#Sum abundances by compartment and phylum
taxa_counterHost <- function(taxa, data, pattern1, pattern2, pattern3) {
  #Pattern 1 and 2 are strings for regex that defines stuff that is unknown.
  outdf <- list()
  k <- 1
  names(data)[grep(pattern1, names(data))] <- "Unknown"
  names(data)[grep(pattern2, names(data))] <- "Unknown"
  names(data)[grep(pattern3, names(data))] <- "Unknown"
  
  for (i in unique(taxa)) {
    if (length(grep(i, names(data))) > 1) {
      df <- aggregate(rowSums(data[, grep(i, names(data))]) ~
                        data$taxon_final + data$compartment, FUN = median)
      df$taxon <- i
      names(df) <- c("host","compartment", "sum", "taxon")
      outdf[[k]] <- df
    } else if (length(grep(i, names(data))) == 1) {
      df <- aggregate(data[, grep(i, names(data))] ~
                        data$taxon_final + data$compartment, FUN = median)
      df$taxon <- i
      names(df) <- c("host","compartment", "sum", "taxon")
      outdf[[k]] <- df
    }
    k <- k + 1
  }
  return(outdf)
}


phyla_counts <- taxa_counterHost(taxa = phyla, 
                             data = merged_dat, 
                             pattern1 = ".*d:Fungi$",
                             pattern2 = "Zotu\\d+_$",
                             pattern3 = "Zotu\\d+_$")

df <- ldply (phyla_counts, data.frame)
df$group <- paste(df$host, df$compartment, sep = "")

library(data.table)
setDT(df)   # coerce to data.table
data_wide <- dcast(df, group ~ taxon, 
                   value.var = c("sum"))

names(data_wide)[-1] == fornext
data_wide_t <- data.frame(t(data_wide[,2:length(data_wide)]))
names(data_wide_t) <- data_wide$group


pdf(width = 12,height = 9, file = "./visuals/phyla_its_taxon_final_barplot_count.pdf")
# Get the stacked barplot
par(mar = c(3,6,1,1), oma = c(16,3,1,1))
a <- barplot(as.matrix(data_wide_t), 
        col = colors,
        border="white", 
        cex.lab = 1.5,
        cex.axis = 1.5,
        las = 2, 
        ylim = c(0,3.5),
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

bact_otus <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                     stringsAsFactors = F)

bact_otus[,3:length(bact_otus)] <- bact_otus[,3:length(bact_otus)] - 1
bact_otus[,3:length(bact_otus)] <- bact_otus[,3:length(bact_otus)] / bact_otus$ISD

bact_meta <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                     stringsAsFactors = F)

bact_otus$sample <- gsub("X", "", bact_otus$sample)

merged_dat <- merge(bact_meta, bact_otus, by.x = "sample", by.y = "sample", all.y = T)

#recode the names of merged dat for each Zotu as the taxonomy hypothesis
for(i in 204:length(merged_dat)){
  if(any(grepl(names(merged_dat)[i], bacts$V1))){
    names(merged_dat)[i] <- paste(names(merged_dat)[i],
                                  bacts$V4[bacts$V1 == names(merged_dat)[i]],
                                  sep = "_")
  }
}

#Extract unique phyla
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

#Sum abundances by compartment and phylum
taxa_counterHost <- function(taxa, data, pattern1, pattern2, pattern3) {
  #Pattern 1 and 2 are strings for regex that defines stuff that is unknown.
  outdf <- list()
  k <- 1
  names(data)[grep(pattern1, names(data))] <- "Unknown"
  names(data)[grep(pattern2, names(data))] <- "Unknown"
  names(data)[grep(pattern3, names(data))] <- "Unknown"
  
  for (i in unique(taxa)) {
    if (length(grep(i, names(data))) > 1) {
      df <- aggregate(rowSums(data[, grep(i, names(data))]) ~
                        data$taxon_final + data$compartment, FUN = median)
      df$taxon <- i
      names(df) <- c("host","compartment", "sum", "taxon")
      outdf[[k]] <- df
    } else if (length(grep(i, names(data))) == 1) {
      df <- aggregate(data[, grep(i, names(data))] ~
                        data$taxon_final + data$compartment, FUN = median)
      df$taxon <- i
      names(df) <- c("host","compartment", "sum", "taxon")
      outdf[[k]] <- df
    }
    k <- k + 1
  }
  return(outdf)
}

phyla_counts <- taxa_counterHost(taxa = phyla, 
                                 data = merged_dat, 
                                 pattern1 = ".*k:Bacteria$",
                                 pattern2 = "Zotu\\d+_$",
                                 pattern3 = "Zotu\\d+_$")

df <- ldply (phyla_counts, 
             data.frame)
df$group <- paste(df$host, 
                  df$compartment, 
                  sep = "")

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

pdf(width = 12,height = 9, file = "./visuals/phyla_16s_taxon_final_barplot_countData.pdf")
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
             ylim = c(0,10),
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
     0, 
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


##########
# legend #
##########

pdf(width = 3,height = 14, file = "./visuals/phyla_16s_barplot_legend_modeledDatatest.pdf")
par(oma = c(10,2,10,1))
plot(x = seq(1,10,1), y = seq(1,10,1),
     type = "n",
     frame.plot = F,
     axes = F,
     xlab = "",
     ylab = "")
legend(legend = rev(rownames(wide)),
       x = "left",
       col = rev(colors_new),
       xpd = NA,
       bty = "n",
       pch = 15)
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
taxa_counter <- function(taxa, data, pattern1, pattern2, pattern3) {
  #Pattern 1 and 2 are strings for regex that defines stuff that is unknown.
  outdf <- list()
  k <- 1
  names(data)[grep(pattern1, names(data))] <- "Unknown"
  names(data)[grep(pattern2, names(data))] <- "Unknown"
  names(data)[grep(pattern3, names(data))] <- "Unknown"
  
  for (i in unique(taxa)) {
    if (length(grep(i, names(data))) > 1) {
      df <- aggregate(rowSums(data[, grep(i, names(data))]) ~
                        data$compartment, FUN = median)
      df$taxon <- i
      names(df) <- c("compartment", "sum", "taxon")
      outdf[[k]] <- df
    } else if (length(grep(i, names(data))) == 1) {
      df <- aggregate(data[, grep(i, names(data))] ~
                        data$compartment, FUN = median)
      df$taxon <- i
      names(df) <- c("compartment", "sum", "taxon")
      outdf[[k]] <- df
    }
    k <- k + 1
  }
  return(outdf)
}

phyla_counts <- taxa_counter(taxa = phyla, 
                                 data = merged_dat, 
                                 pattern1 = ".*k:Bacteria$",
                                 pattern2 = "Zotu\\d+_$",
                                 pattern3 = ".*k:Archaea$")

df <- ldply (phyla_counts, data.frame)

library(tidyverse)
wide = df %>% 
  spread(compartment, sum)
wide
wide <- wide[match( fornext,wide$taxon),]

pdf(width = 5,height = 6, file = "./visuals/phyla_16s_compartment_barplot_countData.pdf")
# Get the stacked barplot
par(mar = c(3,6,1,1), oma = c(2,2,2,2))
barplot(as.matrix(wide[,2:3]), 
        col=colors_new, 
        border="white", 
        cex.lab = 1.5,
        cex.axis = 1.5,
        las = 2, 
        ylim = c(0,0.05),
        space=0.04, 
        font.axis=2, 
        xaxt = "n",
        xlab="")
text( c("EN", "EP"), cex = 1.5, x = c(0.5,1.5), y = -0.005, xpd = NA)
title(ylab = "Abundance (ratio with ISD)", line = 5, cex.lab = 2, xpd = NA)
dev.off()


##########
# legend #
##########
wide$taxon[wide$taxon == "k:Archaea"] <- "Archaea"

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
