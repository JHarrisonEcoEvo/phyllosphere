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

merged_dat16s <- merged_dat16s[-grep("Unkn",merged_dat16s$taxon_final),]

# head(names(merged_dat16s),206)
# tail(names(merged_dat16s),206)
# names(merged_dat16s)[length(merged_dat16s)-4]

#Combine by host first
comboHost <- data.frame(matrix(nrow = 1, ncol = 1+ length(202:(length(merged_dat16s)-4))))

k <- 1
for(i in unique(merged_dat16s$taxon_final)){
  comboHost[k,] <- c(i,colSums(merged_dat16s[merged_dat16s$taxon_final == i,
                                             202:(length(merged_dat16s)-4)]))
  k <- k + 1
}
names(comboHost) <- c("sample", names(merged_dat16s)[ 202:(length(merged_dat16s)-4)])

for(i in 2:length(comboHost)){
  comboHost[,i] <- as.numeric(comboHost[,i])  
}

comboHost$sample <- gsub(" var. .*","", comboHost$sample)

datr16s <-  vegan::decostand(comboHost[,2:length(comboHost)], 
                             method = "hellinger")

dat_e16s <- dist(x  = as.matrix(datr16s), method = "euclidean", upper = T)
labels(dat_e16s)


datits <- read.csv("./processedData/ITSp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv",
                   stringsAsFactors = F)

metadatits <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                       stringsAsFactors = F)

#match order
merged_datits <- merge(metadatits, datits, by.x = "sample", by.y = "sample", all.y = T)
dim(merged_datits)

merged_datits <- merged_datits[-grep("Unkn",merged_datits$taxon_final),]

head(names(merged_datits),206)
tail(names(merged_datits),206)
names(merged_datits)[length(merged_datits)-2]

#Combine by host first
comboHostits <- data.frame(matrix(nrow = 1, ncol = 1+ length(202:(length(merged_datits)-2))))

k <- 1
for(i in unique(merged_datits$taxon_final)){
  comboHostits[k,] <- c(i,colSums(merged_datits[merged_datits$taxon_final == i,
                                                202:(length(merged_datits)-2)]))
  k <- k + 1
}
names(comboHostits) <- c("sample", names(merged_datits)[ 202:(length(merged_datits)-2)])

for(i in 2:length(comboHostits)){
  comboHostits[,i] <- as.numeric(comboHostits[,i])  
}

datrits <-  vegan::decostand(comboHostits[,2:length(comboHostits)], 
                             method = "hellinger")

datrits <- dist(x  = as.matrix(datrits), method = "euclidean", upper = T)

mantel(dat_e16s, datrits)