
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

dendBact <- hclust(dat_e16s, method="average" ) %>% as.dendrogram()

correctOrder <- labels(dendBact)

#fixing error
metadat16s$family.x[metadat16s$taxon_final== "Populus tremulloides"]<- "Salicaceae"
metadat16s$family.x[metadat16s$taxon_final== "Juncus sp."]<- "Juncaceae"

k <- 1
fam <- NA
for(i in comboHost$sample[as.numeric(correctOrder)]){
  print(i)
  if( length(unique(metadat16s$family.x[grep(i, metadat16s$taxon_final)])) > 1){
    print("hey")
  }
  fam[k] <- unique(metadat16s$family.x[grep(i, metadat16s$taxon_final)])
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
plot(dendBact)


##########
# Do over for fungi
rm(list=setdiff(ls(), c("dendBact", "df")))

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
pdf(width = 11, height = 11, file = "./visuals/tanglegram.pdf")
tanglegram(dlist,
           lab.cex = 1,
           common_subtrees_color_lines = FALSE, highlight_distinct_edges  = TRUE, highlight_branches_lwd=FALSE, 
           margin_inner=12,
           lwd=2
)
dev.off()
