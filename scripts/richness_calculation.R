###############################
# Make rarified abundance data #
rm(list=ls())

dat <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG", stringsAsFactors = F)
dat[,3:length(dat)]  <- dat[,3:length(dat)] -1
#dat[,3:length(dat)]  <- dat[,3:length(dat)] / dat$ISD
names(dat)[length(dat)-4]

#devtools::install_github("adw96/breakaway")
library(breakaway)

newmat <- data.frame(t(dat[,3:(length(dat)-4)]))
names(newmat) <- dat$sample
row.names(newmat) <- names(sample)[3:(length(dat)-4)]
newmat <- newmat[rowSums(newmat) != 0,colSums(newmat) != 0]

frequencytablelist <- build_frequency_count_tables(newmat)

rich <- NA
for(i in 1:length(frequencytablelist)){
  richest <- breakaway(frequencytablelist[[i]])
  rich[i] <- richest$estimate
}

names(rich) <- gsub("X","",names(newmat))

meta <- read.csv("processedData/16smetadat_wrangled_for_post_modeling_analysis.csv")
meta <- meta[match(names(rich), meta$sample),]

table(names(rich) == meta$sample)

a <- aov(log(rich) ~ paste(meta$habit,meta$compartment))
xtable::xtable(TukeyHSD(a)$`paste(meta$habit, meta$compartment)`)
# % latex table generated in R 4.1.2 by xtable 1.8-4 package
# % Fri Nov 12 15:30:15 2021
# \begin{table}[ht]
# \centering
# \begin{tabular}{rrrrr}
# \hline
# & diff & lwr & upr & p adj \\ 
# \hline
# forb EP-forb EN & 0.76 & 0.37 & 1.15 & 0.00 \\ 
# graminoid EN-forb EN & -0.46 & -1.31 & 0.40 & 0.74 \\ 
# graminoid EP-forb EN & 0.74 & -0.00 & 1.48 & 0.05 \\ 
# shrub EN-forb EN & 0.14 & -0.36 & 0.63 & 0.99 \\ 
# shrub EP-forb EN & 0.82 & 0.37 & 1.28 & 0.00 \\ 
# tree EN-forb EN & 0.18 & -0.31 & 0.67 & 0.96 \\ 
# tree EP-forb EN & 0.74 & 0.26 & 1.23 & 0.00 \\ 
# graminoid EN-forb EP & -1.22 & -2.06 & -0.38 & 0.00 \\ 
# graminoid EP-forb EP & -0.02 & -0.75 & 0.71 & 1.00 \\ 
# shrub EN-forb EP & -0.62 & -1.10 & -0.15 & 0.00 \\ 
# shrub EP-forb EP & 0.07 & -0.36 & 0.49 & 1.00 \\ 
# tree EN-forb EP & -0.58 & -1.05 & -0.12 & 0.00 \\ 
# tree EP-forb EP & -0.02 & -0.48 & 0.44 & 1.00 \\ 
# graminoid EP-graminoid EN & 1.20 & 0.15 & 2.25 & 0.01 \\ 
# shrub EN-graminoid EN & 0.59 & -0.30 & 1.48 & 0.48 \\ 
# shrub EP-graminoid EN & 1.28 & 0.41 & 2.15 & 0.00 \\ 
# tree EN-graminoid EN & 0.63 & -0.26 & 1.52 & 0.38 \\ 
# tree EP-graminoid EN & 1.20 & 0.31 & 2.09 & 0.00 \\ 
# shrub EN-graminoid EP & -0.61 & -1.40 & 0.18 & 0.28 \\ 
# shrub EP-graminoid EP & 0.08 & -0.68 & 0.85 & 1.00 \\ 
# tree EN-graminoid EP & -0.56 & -1.35 & 0.22 & 0.37 \\ 
# tree EP-graminoid EP & 0.00 & -0.78 & 0.79 & 1.00 \\ 
# shrub EP-shrub EN & 0.69 & 0.16 & 1.22 & 0.00 \\ 
# tree EN-shrub EN & 0.04 & -0.52 & 0.60 & 1.00 \\ 
# tree EP-shrub EN & 0.61 & 0.05 & 1.16 & 0.02 \\ 
# tree EN-shrub EP & -0.65 & -1.17 & -0.12 & 0.00 \\ 
# tree EP-shrub EP & -0.08 & -0.60 & 0.44 & 1.00 \\ 
# tree EP-tree EN & 0.57 & 0.01 & 1.12 & 0.04 \\ 
# \hline
# \end{tabular}
# \end{table}

pdf(file = "visuals/richness-16s.pdf", width = 8, height = 6)
par(oma = c(5,1,1,1))

boxplot((rich) ~ meta$habit + meta$compartment, las = 2, ylab = "estimated richness", xlab = "", main = "16S", outline = F)


text("16S", x = 0, y = 250, xpd = NA, cex= 2)
dev.off()
##########################
rm(list=ls())
library(breakaway)
dat <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG", stringsAsFactors = F)
dat[,3:length(dat)]  <- dat[,3:length(dat)] -1
#dat[,3:length(dat)]  <- dat[,3:length(dat)] / dat$ISD
names(dat)[length(dat)-2]

newmat <- data.frame(t(dat[,3:(length(dat)-4)]))
names(newmat) <- dat$sample
row.names(newmat) <- names(sample)[3:(length(dat)-4)]
newmat <- newmat[rowSums(newmat) != 0,colSums(newmat) != 0]

frequencytablelist <- build_frequency_count_tables(newmat)

rich <- NA
for(i in 1:length(frequencytablelist)){
  richest <- breakaway(frequencytablelist[[i]])
  rich[i] <- richest$estimate
}

names(rich) <- gsub("X","",names(newmat))

meta <- read.csv("processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv")
meta <- meta[match(names(rich), meta$sample),]

table(names(rich) == meta$sample)

a <- aov(log(rich) ~ paste(meta$habit,meta$compartment))
TukeyHSD(a)
xtable::xtable(TukeyHSD(a)$`paste(meta$habit, meta$compartment)`)

pdf(file = "visuals/richness-its.pdf", width = 8, height = 6)
par(oma = c(5,1,1,1))
boxplot((rich) ~ meta$habit + meta$compartment, las = 2, ylab = "estimated richness", xlab = "", main = "", outline = F)

text("ITS", x = 0, y = 250, xpd = NA, cex= 2)
dev.off()
