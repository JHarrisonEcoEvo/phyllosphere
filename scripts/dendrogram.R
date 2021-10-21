
###########################
# cluster dendrogram#
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

dat16s <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                   stringsAsFactors = F)

dat16s$sample <- gsub("X","",dat16s$sample)
dat16s[,grep("Zotu", names(dat16s))] <- dat16s[,grep("Zotu", names(dat16s))] - 1
dat16s[,grep("Zotu", names(dat16s))] <- dat16s[,grep("Zotu", names(dat16s))] / dat16s$ISD

metadat16s <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                       stringsAsFactors = F)

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

#Subset to a specific compartment
merged_dat16s_en <- merged_dat16s[merged_dat16s$compartment == "EP",]

#Combine by host first
comboHost <- data.frame(matrix(nrow = 1, ncol = 1+ length(205:(length(merged_dat16s_en)-4))))

merged_dat16s_en$taxon_final <- gsub(" var. .*","", merged_dat16s_en$taxon_final)

k <- 1
for(i in unique(merged_dat16s_en$taxon_final)){
  comboHost[k,] <- c(i,colSums(merged_dat16s_en[merged_dat16s_en$taxon_final == i,
                                 205:(length(merged_dat16s_en)-4)]))
  k <- k + 1
}
names(comboHost) <- c("sample", names(merged_dat16s_en)[ 205:(length(merged_dat16s_en)-4)])

for(i in 2:length(comboHost)){
  comboHost[,i] <- as.numeric(comboHost[,i])  
}


datr16s <-  vegan::decostand(comboHost[,2:length(comboHost)], 
                             method = "hellinger")

dat_e16s <- dist(x  = as.matrix(datr16s), method = "euclidean", upper = T)
labels(dat_e16s)

dendBact <- hclust(dat_e16s, method="average" ) %>% as.dendrogram()

correctOrder <- labels(dendBact)

#fixing error
merged_dat16s_en$family.x[merged_dat16s_en$taxon_final== "Populus tremulloides"]<- "Salicaceae"
merged_dat16s_en$family.x[merged_dat16s_en$taxon_final== "Juncus sp."]<- "Juncaceae"

k <- 1
fam <- NA
for(i in comboHost$sample[as.numeric(correctOrder)]){
  print(i)
  if( length(unique(merged_dat16s_en$family.x[grep(i, merged_dat16s_en$taxon_final)])) > 1){
    print("hey")
  }
  fam[k] <- unique(merged_dat16s_en$family.x[grep(i, merged_dat16s_en$taxon_final)])
  k<- k + 1
}
#sanity check
cbind(comboHost$sample[as.numeric(correctOrder)],
      fam)
#assign colors by family
skittles <- rainbow(n = length(unique(fam)))

colorsHost <- rep("red", length(fam))

k <- 1
for(i in unique(fam)){
  colorsHost[grep(i, fam)] <- skittles[k]
  k <- k + 1
}


#sanity check
df <- data.frame(comboHost$sample[as.numeric(correctOrder)],
      fam, colorsHost)


#Plot a legend
pdf(width = 4, height = 12, file = "./visuals/tanglegram_legend.pdf")
plot(NULL,ylim = c(0,200),
     xlim = c(1,10),
     xaxt = "n",
     yaxt = "n",
     frame.plot = F, 
     ylab  = "",
     xlab = ""
     )
  legend(x = "topleft",
         legend = as.character(unique(df$fam)),
       col = as.character(unique(df$colorsHost)),
       cex = 1,
       pch = 16)
dev.off()

dendBact <- dendBact %>% 
  # Custom branches
  set("branches_col", "grey") %>% set("branches_lwd", 3) %>%
  # Custom labels
  set("labels", comboHost$sample[as.numeric(correctOrder)]) %>% set("labels_cex", 0.8) %>%
  set("labels_colors", colorsHost)
#testplot for colors to make sure they match the tanglgram 
#plot(dendBact)


##########
# Do over for fungi
rm(list=setdiff(ls(), c("dendBact", "df", "metadatits")))

datits <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                   stringsAsFactors = F)

datits$sample <- gsub("X","",datits$sample)
datits[,grep("Zotu", names(datits))] <- datits[,grep("Zotu", names(datits))] - 1
datits[,grep("Zotu", names(datits))] <- datits[,grep("Zotu", names(datits))] / datits$ISD

#match order
merged_datits <- merge(metadatits, datits, by.x = "sample", by.y = "sample")
dim(merged_datits)

#Subset to a specific compartment
merged_datits <- merged_datits[merged_datits$compartment == "EP",]

#merged_datits <- merged_datits[-grep("Unkn",merged_datits$taxon_final),]

# head(names(merged_datits),206)
# tail(names(merged_datits),206)
# names(merged_datits)[length(merged_datits)-2]

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

datrits <-  vegan::decostand(comboHostits[,2:length(comboHostits)], 
                             method = "hellinger")

datrits <- dist(x  = as.matrix(datrits), method = "euclidean", upper = T)
labels(datrits)

dendFung <- hclust(datrits, method="average" ) %>% as.dendrogram()

correctOrderF <- labels(dendFung)

#remove varieties for cleaner ploting
comboHostits$sample <- gsub(" var. .*","", comboHostits$sample)

#fixing error
metadatits$family.x[metadatits$taxon_final== "Populus tremulloides"]<- "Salicaceae"
metadatits$family.x[metadatits$taxon_final== "Juncus sp."]<- "Juncaceae"


k <- 1
fam <- NA
for(i in comboHostits$sample[as.numeric(correctOrderF)]){
  print(i)
  if( length(unique(comboHostits$family.x[grep(i, comboHostits$taxon_final)])) > 1){
    print("hey")
  }
  fam[k] <- unique(metadatits$family.x[grep(i, metadatits$taxon_final)])
  k<- k + 1
}
#sanity check
cbind(comboHostits$sample[as.numeric(correctOrderF)],
      fam)
#assign colors by family
#We use the same colormatches as we did before between family and color. 
#this convolution is to ensure we have the same matches
colorsHost <- rep("red", length(fam))

df2 <- data.frame(comboHostits$sample[as.numeric(correctOrderF)],
             fam)
df2$color <- "red" #placeholder

for(i in 1:length(df2$fam)){ #loop by item in family in the new dataframe
  #find the color in the old dataframe thatmatches the new family and stick it into df2
  df2$color[i] <- as.character(unique(df$colorsHost[df$fam == df2$fam[i]]))

  k <- k + 1
}

#sanity check
cbind(comboHostits$sample[as.numeric(correctOrderF)],
      df2)

dendFung <- dendFung %>% 
  # Custom branches
  set("branches_col", "grey") %>% set("branches_lwd", 3) %>%
  # Custom labels
  set("labels", comboHostits$sample[as.numeric(correctOrderF)]) %>% set("labels_cex", 0.8) %>%
  set("labels_colors", df2$color)

#par(mar=c(14, 3,2,2), oma = c(4,1,1,1))
plot(dendFung, xpd =NA, yaxt = "n")

######
#try tanglegram https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html#comparing-two-dendrograms
#dendlist chains the dendrograms together to make them easier to work with
dlist <-dendlist(dendBact, dendFung)
# 
# untangled <- dlist %>% untangle( method = "step2side")
# 
# tanglegram(untangled,
#            common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
#            margin_inner=7,
#            lwd=2
# )
#pdf(width = 11, height = 11, file = "./visuals/tanglegram_en.pdf")
pdf(width = 11, height = 11, file = "./visuals/tanglegram_ep.pdf")

tanglegram(dlist,
           lab.cex = 1,
           common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
           margin_inner=12,
           lwd=2
)
dev.off()

# 
# ##########################################
# # Within taxon comparison by compartment #
# ##########################################
# 
# rm(list = ls())
# library(ecodist)
# library(vegan)
# library(dendextend)
# 
# #Function from Mage
# add.alpha <- function(col, alpha=1){
#   if(missing(col))
#     stop("Please provide a vector of colours.")
#   apply(sapply(col, col2rgb)/255, 2, 
#         function(x) 
#           rgb(x[1], x[2], x[3], alpha=alpha))  
# }
# 
# dat16s <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
#                    stringsAsFactors = F)
# 
# dat16s$sample <- gsub("X","",dat16s$sample)
# dat16s[,grep("Zotu", names(dat16s))] <- dat16s[,grep("Zotu", names(dat16s))] - 1
# dat16s[,grep("Zotu", names(dat16s))] <- dat16s[,grep("Zotu", names(dat16s))] / dat16s$ISD
# 
# metadat16s <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
#                        stringsAsFactors = F)
# 
# metadatits <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
#                        stringsAsFactors = F)
# metadat16s <-
#   metadat16s[metadat16s$sample%in%  metadatits$sample,] 
# table(metadat16s$sample  %in% metadatits$sample)
# 
# metadatits <-
#   metadatits[metadatits$sample %in% metadat16s$sample,] 
# table(metadatits$sample %in% metadat16s$sample)
# 
# 
# #match order
# merged_dat16s <- merge(metadat16s, dat16s, by.x = "sample", by.y = "sample")
# dim(merged_dat16s)
# 
# # head(names(merged_dat16s),206)
# # tail(names(merged_dat16s),206)
# # names(merged_dat16s)[length(merged_dat16s)-4]
# 
# #Subset to a specific compartment
# merged_dat16s_en <- merged_dat16s[merged_dat16s$compartment == "EP",]
# 
# #Combine by host first
# comboHost <- data.frame(matrix(nrow = 1, ncol = 1+ length(205:(length(merged_dat16s_en)-4))))
# 
# merged_dat16s_en$taxon_final <- gsub(" var. .*","", merged_dat16s_en$taxon_final)
# 
# k <- 1
# for(i in unique(merged_dat16s_en$taxon_final)){
#   comboHost[k,] <- c(i,colSums(merged_dat16s_en[merged_dat16s_en$taxon_final == i,
#                                                 205:(length(merged_dat16s_en)-4)]))
#   k <- k + 1
# }
# names(comboHost) <- c("sample", names(merged_dat16s_en)[ 205:(length(merged_dat16s_en)-4)])
# 
# for(i in 2:length(comboHost)){
#   comboHost[,i] <- as.numeric(comboHost[,i])  
# }
# 
# 
# datr16s <-  vegan::decostand(comboHost[,2:length(comboHost)], 
#                              method = "hellinger")
# 
# dat_e16s <- dist(x  = as.matrix(datr16s), method = "euclidean", upper = T)
# labels(dat_e16s)
# 
# dendBact <- hclust(dat_e16s, method="average" ) %>% as.dendrogram()
# 
# correctOrder <- labels(dendBact)
# 
# #fixing error
# merged_dat16s_en$family.x[merged_dat16s_en$taxon_final== "Populus tremulloides"]<- "Salicaceae"
# merged_dat16s_en$family.x[merged_dat16s_en$taxon_final== "Juncus sp."]<- "Juncaceae"
# 
# k <- 1
# fam <- NA
# for(i in comboHost$sample[as.numeric(correctOrder)]){
#   print(i)
#   if( length(unique(merged_dat16s_en$family.x[grep(i, merged_dat16s_en$taxon_final)])) > 1){
#     print("hey")
#   }
#   fam[k] <- unique(merged_dat16s_en$family.x[grep(i, merged_dat16s_en$taxon_final)])
#   k<- k + 1
# }
# #sanity check
# cbind(comboHost$sample[as.numeric(correctOrder)],
#       fam)
# #assign colors by family
# skittles <- rainbow(n = length(unique(fam)))
# 
# colorsHost <- rep("red", length(fam))
# 
# k <- 1
# for(i in unique(fam)){
#   colorsHost[grep(i, fam)] <- skittles[k]
#   k <- k + 1
# }
# 
# 
# #sanity check
# df <- data.frame(comboHost$sample[as.numeric(correctOrder)],
#                  fam, colorsHost)
# 
# 
# dendBactep <- dendBact %>% 
#   # Custom branches
#   set("branches_col", "grey") %>% set("branches_lwd", 3) %>%
#   # Custom labels
#   set("labels", comboHost$sample[as.numeric(correctOrder)]) %>% set("labels_cex", 0.8) %>%
#   set("labels_colors", colorsHost)
# #testplot for colors to make sure they match the tanglgram 
# #plot(dendBact)
# 
# #EN
# 
# dat16s <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
#                    stringsAsFactors = F)
# 
# dat16s$sample <- gsub("X","",dat16s$sample)
# dat16s[,grep("Zotu", names(dat16s))] <- dat16s[,grep("Zotu", names(dat16s))] - 1
# dat16s[,grep("Zotu", names(dat16s))] <- dat16s[,grep("Zotu", names(dat16s))] / dat16s$ISD
# 
# metadat16s <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
#                        stringsAsFactors = F)
# 
# metadatits <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
#                        stringsAsFactors = F)
# metadat16s <-
#   metadat16s[metadat16s$sample%in%  metadatits$sample,] 
# table(metadat16s$sample  %in% metadatits$sample)
# 
# metadatits <-
#   metadatits[metadatits$sample %in% metadat16s$sample,] 
# table(metadatits$sample %in% metadat16s$sample)
# 
# 
# #match order
# merged_dat16s <- merge(metadat16s, dat16s, by.x = "sample", by.y = "sample")
# dim(merged_dat16s)
# 
# # head(names(merged_dat16s),206)
# # tail(names(merged_dat16s),206)
# # names(merged_dat16s)[length(merged_dat16s)-4]
# 
# #Subset to a specific compartment
# merged_dat16s_en <- merged_dat16s[merged_dat16s$compartment == "EN",]
# 
# #Combine by host first
# comboHost <- data.frame(matrix(nrow = 1, ncol = 1+ length(205:(length(merged_dat16s_en)-4))))
# 
# merged_dat16s_en$taxon_final <- gsub(" var. .*","", merged_dat16s_en$taxon_final)
# 
# k <- 1
# for(i in unique(merged_dat16s_en$taxon_final)){
#   comboHost[k,] <- c(i,colSums(merged_dat16s_en[merged_dat16s_en$taxon_final == i,
#                                                 205:(length(merged_dat16s_en)-4)]))
#   k <- k + 1
# }
# names(comboHost) <- c("sample", names(merged_dat16s_en)[ 205:(length(merged_dat16s_en)-4)])
# 
# for(i in 2:length(comboHost)){
#   comboHost[,i] <- as.numeric(comboHost[,i])  
# }
# 
# 
# datr16s <-  vegan::decostand(comboHost[,2:length(comboHost)], 
#                              method = "hellinger")
# 
# dat_e16s <- dist(x  = as.matrix(datr16s), method = "euclidean", upper = T)
# labels(dat_e16s)
# 
# dendBact <- hclust(dat_e16s, method="average" ) %>% as.dendrogram()
# 
# correctOrder <- labels(dendBact)
# 
# #fixing error
# merged_dat16s_en$family.x[merged_dat16s_en$taxon_final== "Populus tremulloides"]<- "Salicaceae"
# merged_dat16s_en$family.x[merged_dat16s_en$taxon_final== "Juncus sp."]<- "Juncaceae"
# 
# k <- 1
# fam <- NA
# for(i in comboHost$sample[as.numeric(correctOrder)]){
#   print(i)
#   if( length(unique(merged_dat16s_en$family.x[grep(i, merged_dat16s_en$taxon_final)])) > 1){
#     print("hey")
#   }
#   fam[k] <- unique(merged_dat16s_en$family.x[grep(i, merged_dat16s_en$taxon_final)])
#   k<- k + 1
# }
# #sanity check
# cbind(comboHost$sample[as.numeric(correctOrder)],
#       fam)
# #assign colors by family
# #We use the same colormatches as we did before between family and color. 
# #this convolution is to ensure we have the same matches
# colorsHost <- rep("red", length(fam))
# 
# df2 <- data.frame(comboHostits$sample[as.numeric(correctOrderF)],
#                   fam)
# df2$color <- "red" #placeholder
# 
# for(i in 1:length(df2$fam)){ #loop by item in family in the new dataframe
#   #find the color in the old dataframe thatmatches the new family and stick it into df2
#   df2$color[i] <- as.character(unique(df$colorsHost[df$fam == df2$fam[i]]))
#   
#   k <- k + 1
# }
# 
# #sanity check
# cbind(comboHostits$sample[as.numeric(correctOrderF)],
#       df2)
# 
# #sanity check
# df <- data.frame(comboHost$sample[as.numeric(correctOrder)],
#                  fam, colorsHost)
# 
# 
# dendBacten <- dendBact %>% 
#   # Custom branches
#   set("branches_col", "grey") %>% set("branches_lwd", 3) %>%
#   # Custom labels
#   set("labels", comboHost$sample[as.numeric(correctOrder)]) %>% set("labels_cex", 0.8) %>%
#   set("labels_colors", colorsHost)
# #testplot for colors to make sure they match the tanglgram 
# #plot(dendBact)
# 
# dlist <-dendlist(dendBactep, dendBacten)
#  
# # untangled <- dlist %>% untangle( method = "step2side")
# # 
# # tanglegram(untangled,
# #            common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
# #            margin_inner=7,
# #            lwd=2
# # )
# 
# pdf(width = 11, height = 11, file = "./visuals/tanglegram_bacteria_enVsep.pdf")
# 
# tanglegram(dlist,
#            lab.cex = 1,
#            common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
#            margin_inner=12,
#            lwd=2
# )
# dev.off()

