rm(list=ls())

#Bring in results and figure out which taxa were predicted well enough to warrant
#extraction of important features. 

results <- read.csv("modelingResults/all_its_occupancy.csv", stringsAsFactors = F)
results <- results[results$taxon != "taxon",]
results <- results[as.numeric(results$mcc_nested_resampling) > 0,]

#Bring in variable importance metrics for only those taxa that were predicted
dat <- list.files("modelingResults/varImp_occupancy/")
dat <- dat[grep("*ITS*", dat)]
#dat <- dat[-grep("*Reduced*", dat)]
#dat <- dat[-grep("*NOHOST*", dat)]

dat <- grep(paste(results$taxon, collapse="|"), 
             dat, value=TRUE)

varimps <- lapply(paste("modelingResults/varImp_occupancy/", dat, sep = ""), read.csv)
length(varimps)

#Extract top ten from each and then do table to see which are most useful. 

#look at this clunker!

important <- vector()

for(i in 1:length(varimps)){
  important <- c(important,
               as.character(varimps[[i]][head(rev(order(varimps[[i]][,2])), n = 10),1]))
}
table(important)
sort(table(important))
hist(sort(table(important)))
length(varimps)

#Bring in taxonomy for each taxon and see if there are certain patterns by phyla
#where, say, a certain suite of predictors is more important for certain phyla
#If this appears to be the case, then make a heatmap figure


library(stringr)

#Function by F. Bredoire
read.sintax <- function(x) {
  d <- read.table(x, sep = "\t", header = FALSE, fill = TRUE, na.strings = "")
  d$OTU_ID <- str_extract(d$V1, "otu[0-9]+")
  d$DOMAIN <- str_extract(d$V4, "d:[^,]+")
  d$KINGDOM <- str_extract(d$V4, "k:[^,]+")
  d$PHYLUM <- str_extract(d$V4, "p:[^,]+")
  d$CLASS <- str_extract(d$V4, "c:[^,]+")
  d$ORDER <- str_extract(d$V4, "o:[^,]+")
  d$FAMILY <- str_extract(d$V4, "f:[^,]+")
  d$GENUS <- str_extract(d$V4, "g:[^,]+")
  d$SPECIES <- str_extract(d$V4, "s:[^,]+")
  return(d[, 5:13])
}

tax <- read.sintax("processedData/smallmem97_ITS.sintax")
tax$OTU_ID <- gsub("otu", "Zotu", tax$OTU_ID)

#subset to just those taxa that were modeled. 
#calculate table for important features for each taxon. 
#figure out how to display, if interesting.

hits <- gsub("variableImportan.*(Zotu\\d+ITS).*", "\\1", dat)
hits <- gsub("ITS*", "", hits)
tax <- tax[tax$OTU_ID %in% hits,]

names(varimps) <- hits

#ascos
ascos <- na.omit(tax$OTU_ID[tax$PHYLUM == "p:Ascomycota"])

queryindices <- which(names(varimps) %in% ascos)

importantAscos <- vector()
for(i in queryindices){
  importantAscos <- c(important,
                 as.character(varimps[[i]][head(rev(order(varimps[[i]][,2])), n = 10),1]))
}
sort(table(importantAscos))

#basidios
basidios <- na.omit(tax$OTU_ID[tax$PHYLUM == "p:Basidiomycota"])
queryindices <- which(names(varimps) %in% basidios)

importantbasidios <- vector()
for(i in queryindices){
  importantbasidios <- c(important,
                      as.character(varimps[[i]][head(rev(order(varimps[[i]][,2])), n = 10),1]))
}
sort(table(importantbasidios))

#do by class
#subset to just those classes for which we have 5 or more otus


classinfo <- list()
k <- 1
for(j in names(table(tax$CLASS)[table(tax$CLASS) > 5])){
  taxC <- tax[tax$CLASS == j, ]
  taxC <- taxC[!is.na(taxC$CLASS),]
  
  #get the indices for things in that class
  queryindices <- which(names(varimps) %in% taxC$OTU_ID)

  #get top ten most important features
  importantC <- vector()
  for(i in queryindices){
    importantC <- c(importantC,
                         as.character(varimps[[i]][head(rev(order(varimps[[i]][,2])), n = 10),1]))
  }
  #subset to just those features that are important for 25% or more of samples.
  #Those weeding out some features that are only important for a few taxa

  classinfo[[k]] <-  table(importantC)[ table(importantC) > round(0.20*length(queryindices))]
  names(classinfo)[[k]] <- j
  k <- k + 1
}


length(unique(unlist(classinfo)))

# #Need to get the data into something like this, to use the heatmap
#ug. what a pain.
# data <- as.matrix(mtcars)
# 
# # Default Heatmap
# heatmap(data)


nms <- vector()
for(i in 1:length(classinfo)){
  nms <- c(nms, names(classinfo[[i]]))
}

outmat <- list()
k <- 1
for(i in unique(nms)){
  outmat <- c(outmat, i=(lapply(classinfo, FUN=function(x){x[names(x)==i]})))
 # names(outmat[[k]]) <- i
  k <- k + 1
}

#catching and recoding the zero length vectors as zeros
#outmat <- outmat[-c(1,2,3)]

for(i in 1:length(outmat)){
  if(length(outmat[[i]]) == 0){
    outmat[[i]] <- 0
  }
}

forheat <- data.frame(
  matrix(unlist(outmat), nrow = 3 , ncol = length(unique(nms))))
names(forheat) <- unique(nms)
forheat

#REmove host for visualization. We will state clearly in the text that host is important
#thus including it here is redundant
#forheat <- forheat[, -grep("taxon",names(forheat))]

#found this tidbit here, but doesn't work in this case. Saving for posterity:
#https://stackoverflow.com/questions/35089898/transforming-a-list-with-differing-number-of-elements-to-a-data-frame
#forheat <- t(sapply(classinfo, `length<-`, max(lengths(classinfo))))


pdf(width = 10, height = 10, file = "./visuals/its_heatmap_feature_importance_occupancy.pdf")
par(mar=c(3,3,3,3), oma = c(8,3,3,12))
heatmap(as.matrix(forheat),scale = "none",
        col= rev(heat.colors(5)),
        #wesanderson::wes_palette("FantasticFox1", n = 5, type = "discrete"),
        Colv = NA, 
        labCol = c("Leaf area", 
                   "B",
                   "Dead down",
                   "Elevation",
                   "G",
                   "Height sample",
                   "Julian date",
                   "Latitude",
                   "Light intensity (PAR)",
                   
                   "Longitude",
                   "Temp. April",
                   "MEM #2",
                   "Plant volume",
                   "Precip. April",
                   "Pressure",
                   "R",
                   "Shannon's flora",
                   "Shrub richness",
                   "Leaf density (SLA)",
                   "Fm Prime",
                   "Num. leaves extracted",
                   "LEF",
                   "Rel. chlorophyll",
                   "SPAD 420",
                   "Fo prime"
                   ),
        labRow = gsub("c:", "", names(classinfo))
)
dev.off()

pdf(width = 5, height = 5, file = "./visuals/its_heatmap_legend_occupancy.pdf")
plot(NULL)
legend("center",
       bty = "n",
       xpd = NA,
       legend=c("0-4", "5-8", "9-16", "17-22", "22+"),
       fill=c(rev(heat.colors(5))))
#wesanderson::wes_palette("FantasticFox1", n = 5, type = "discrete"))

dev.off()
#ORIGINAL VERSION
#Useful for results of models of data that were not normalized by the ISD. 

# pdf(width = 10, height = 10, file = "./visuals/its_heatmap_feature_importance.pdf")
# par(mar=c(3,3,3,3), oma = c(8,3,3,8))
# heatmap(forheat,scale = "none",
#         col= rev(heat.colors(5)),
#           #wesanderson::wes_palette("FantasticFox1", n = 5, type = "discrete"),
#         Colv = NA, 
#         labCol = c("Compartment EN", 
#                    "Elevation",
#                    "Latitude",
#                    "Num. leaves extracted",
#                    "Mean temp. April",
#                    "Currently flowering",
#                    "Atmospheric pressure",
#                    "Relative chlorophyll",
#                    "Shrub richness",
#                    "leaf density (SLA)",
#                    "SPAD 420",
#                    "Antennaria media (host)",
#                    "Juniperus communis (host)",
#                    "Tree richness"),
#         labRow = gsub("c:", "", row.names(forheat))
#         )
# dev.off()
# 
# pdf(width = 5, height = 5, file = "./visuals/its_heatmap_legend.pdf")
# plot(NULL)
# legend("center",
#        bty = "n",
#        xpd = NA,
#        legend=c("NA", "1-6", "7-9", "10-14", "15", "16+"),
#        fill=c("white",rev(heat.colors(5))))
#          #wesanderson::wes_palette("FantasticFox1", n = 5, type = "discrete"))
# 
# dev.off()

###########################
# Do for bacteria #####
#NOT ENOUGH RESULTS FOR THIS TO BE WORTHWHILE
rm(list=ls())

#Bring in results and figure out which taxa were predicted well enough to warrant
#extraction of important features. 

results <- read.csv("modelingResults/all16S_landscape_count.csv", stringsAsFactors = F)
results <- results[results$rsq_nested_resampling > 0.01,]
results <- results[results$taxon != "taxon",]

#Bring in variable importance metrics for only those taxa that were predicted
dat <- list.files("modelingResults/varImp_landscapeCOUNT/")
dat <- dat[grep("*16S*", dat)]
dat <- dat[-grep("*Reduced*", dat)]

dat <- grep(paste(results$taxon, collapse="|"), 
            dat, value=TRUE)

varimps <- lapply(paste("modelingResults/varImp_landscapeCOUNT/", dat, sep = ""), read.csv)
length(varimps)

#Extract top ten from each and then do table to see which are most useful. 

#look at this clunker!

important <- vector()

for(i in 1:length(varimps)){
  important <- c(important,
                 as.character(varimps[[i]][head(rev(order(varimps[[i]][,2])), n = 10),1]))
}
table(important)
sort(table(important))
hist(sort(table(important)))
length(varimps) #168

#Bring in taxonomy for each taxon and see if there are certain patterns by phyla
#where, say, a certain suite of predictors is more important for certain phyla
#If this appears to be the case, then make a heatmap figure


library(stringr)

#Function by F. Bredoire
read.sintax <- function(x) {
  d <- read.table(x, sep = "\t", header = FALSE, fill = TRUE, na.strings = "")
  d$OTU_ID <- str_extract(d$V1, "otu[0-9]+")
  d$DOMAIN <- str_extract(d$V4, "d:[^,]+")
  d$KINGDOM <- str_extract(d$V4, "k:[^,]+")
  d$PHYLUM <- str_extract(d$V4, "p:[^,]+")
  d$CLASS <- str_extract(d$V4, "c:[^,]+")
  d$ORDER <- str_extract(d$V4, "o:[^,]+")
  d$FAMILY <- str_extract(d$V4, "f:[^,]+")
  d$GENUS <- str_extract(d$V4, "g:[^,]+")
  d$SPECIES <- str_extract(d$V4, "s:[^,]+")
  return(d[, 5:13])
}

tax <- read.sintax("processedData/smallmem97_16S.sintax")
tax$OTU_ID <- gsub("otu", "Zotu", tax$OTU_ID)

#subset to just those taxa that were modeled. 
#calculate table for important features for each taxon. 
#figure out how to display, if interesting.

hits <- gsub("variableImportan.*(Zotu\\d+16S).*", "\\1", dat)
hits <- gsub("16S*", "", hits)
tax <- tax[tax$OTU_ID %in% hits,]

names(varimps) <- hits

#do by class
#subset to just those classes for which we have 5 or more otus


classinfo <- list()
k <- 1
for(j in names(table(tax$CLASS))){
  taxC <- tax[tax$CLASS == j, ]
  taxC <- taxC[!is.na(taxC$CLASS),]
  
  queryindices <- which(names(varimps) %in% taxC$OTU_ID)
  
  importantC <- vector()
  for(i in queryindices){
    importantC <- c(importantC,
                    as.character(varimps[[i]][head(rev(order(varimps[[i]][,2])), n = 10),1]))
  }
  #subset to just those features that are important for 25% or more of samples.
  #Those weeding out some features that are only important for a few taxa
  
  classinfo[[k]] <-  table(importantC)[ table(importantC) > round(0.30*length(queryindices))]
  names(classinfo)[[k]] <- j
  k <- k + 1
}


length(unique(unlist(classinfo)))

# #Need to get the data into something like this, to use the heatmap
# data <- as.matrix(mtcars)
# 
# # Default Heatmap
# heatmap(data)

#found this tidbit here:
#https://stackoverflow.com/questions/35089898/transforming-a-list-with-differing-number-of-elements-to-a-data-frame
forheat <- t(sapply(classinfo, `length<-`, max(lengths(classinfo))))

#pdf(width = 8, height = 8, file = "./visuals/sixteens_heatmap_feature_importance.pdf")
par(mar=c(3,3,3,3), oma = c(8,3,3,8))
heatmap(forheat,scale = "none",
        col= rev(heat.colors(5)),
        #wesanderson::wes_palette("FantasticFox1", n = 5, type = "discrete"),
        Colv = NA, 
        labCol = c("Compartment EN", 
                   "Elevation",
                   "Latitude",
                   "Num. leaves extracted",
                   "Mean temp. April",
                   "Currently flowering",
                   "Atmospheric pressure",
                   "Relative chlorophyll",
                   "Shrub richness",
                   "leaf density (SLA)",
                   "SPAD 420",
                   "Antennaria media (host)",
                   "Juniperus communis (host)",
                   "Tree richness"),
        labRow = gsub("c:", "", row.names(forheat))
)

legend(x=150,
       y = 1,
       xpd = NA,
       legend=c("NA", "1-6", "7-9", "10-14", "15", "16+"),
       fill=c("white",rev(heat.colors(5))))
#wesanderson::wes_palette("FantasticFox1", n = 5, type = "discrete"))

#dev.off()


