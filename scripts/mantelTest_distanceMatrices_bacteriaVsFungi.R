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

dat16s <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                   stringsAsFactors = F)

dat16s$sample <- gsub("X","",dat16s$sample)
dat16s[,grep("Zotu", names(dat16s))] <- dat16s[,grep("Zotu", names(dat16s))] - 1
dat16s[,grep("Zotu", names(dat16s))] <- dat16s[,grep("Zotu", names(dat16s))] / dat16s$ISD

metadat16s <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                       stringsAsFactors = F)
#metadat16s$taxon_final[metadat16s$taxon_final == "Unknown spruce" ] <- "Unknown Spruce"
#write.csv(metadat16s, file = "./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
 #         row.names = F)

metadatits <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                       stringsAsFactors = F)
metadat16s <-
  metadat16s[metadat16s$sample%in%  metadatits$sample,] 
table(metadat16s$sample  %in% metadatits$sample)

metadatits <-
  metadatits[metadatits$sample %in% metadat16s$sample,] 
table(metadatits$sample %in% metadat16s$sample)


#match order
merged_dat16s <- merge(metadat16s, dat16s, by.x = "sample", by.y = "sample")
dim(merged_dat16s)

#merged_dat16s <- merged_dat16s[-grep("Unkn",merged_dat16s$taxon_final),]

# head(names(merged_dat16s),206)
# tail(names(merged_dat16s),206)
# names(merged_dat16s)[length(merged_dat16s)-4]

#Combine by host first
comboHost <- data.frame(matrix(nrow = 1, ncol = 1+ length(205:(length(merged_dat16s)-4))))

k <- 1
for(i in unique(merged_dat16s$taxon_final)){
  comboHost[k,] <- c(i,colSums(merged_dat16s[merged_dat16s$taxon_final == i,
                                             205:(length(merged_dat16s)-4)]))
  k <- k + 1
}
names(comboHost) <- c("sample", names(merged_dat16s)[ 205:(length(merged_dat16s)-4)])

for(i in 2:length(comboHost)){
  comboHost[,i] <- as.numeric(comboHost[,i])  
}

comboHost$sample <- gsub(" var. .*","", comboHost$sample)

datr16s <-  vegan::decostand(comboHost[,2:length(comboHost)], 
                             method = "hellinger")

dat_e16s <- dist(x  = as.matrix(datr16s), method = "euclidean", upper = T)
labels(dat_e16s)
str(dat_e16s)

datits <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                   stringsAsFactors = F)
datits$sample <- gsub("X","",datits$sample)
datits[,grep("Zotu", names(datits))] <- datits[,grep("Zotu", names(datits))] - 1
datits[,grep("Zotu", names(datits))] <- datits[,grep("Zotu", names(datits))] / datits$ISD



#match order
merged_datits <- merge(metadatits, datits, by.x = "sample", by.y = "sample")


#merged_datits <- merged_datits[-grep("Unkn",merged_datits$taxon_final),]

head(names(merged_datits),206)
tail(names(merged_datits),206)
names(merged_datits)[length(merged_datits)-2]

#Combine by host first
comboHostits <- data.frame(matrix(nrow = 1, ncol = 1+ length(205:(length(merged_datits)-2))))

k <- 1
for(i in unique(merged_datits$taxon_final)){
  comboHostits[k,] <- c(i,colSums(merged_datits[merged_datits$taxon_final == i,
                                                205:(length(merged_datits)-2)]))
  k <- k + 1
}
names(comboHostits) <- c("sample", names(merged_datits)[ 205:(length(merged_datits)-2)])

for(i in 2:length(comboHostits)){
  comboHostits[,i] <- as.numeric(comboHostits[,i])  
}

comboHostits$sample <- gsub(" var. .*","", comboHostits$sample)

datrits <-  vegan::decostand(comboHostits[,2:length(comboHostits)], 
                             method = "hellinger")

datrits <- dist(x  = as.matrix(datrits), method = "euclidean", upper = T)
str(datrits)
str(dat_e16s)

mantel(dat_e16s, datrits)