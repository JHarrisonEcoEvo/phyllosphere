#Make plots showing the taxonomy of the taxa that were modeled using random forest
# Groups of taxa to model include:
#1. 16s taxa modeled across all hosts
#2. ITS taxa modeled across all hosts
#1. 16s taxa modeled within hosts, by host. 
#2. ITS taxa modeled within hosts, by host.

#Extras. 
#which taxa were the most prevalent. Make a table. Were they modeled or not. 
#What was distribution of their abundances?

########################################
# 1. 16s taxa modeled across all hosts #
########################################

rm(list=ls())

tax <- read.csv("./processedData/sixteenS_taxa_to_model_via_randomforest.csv")
head(tax)
sintax <- read.table("./processedData/smallmem97_16S.sintax", fill = T)

hits <- sintax[sintax$V1 %in% tax$taxa_16s,]
write.csv(hits, file = "./processedData/taxonomy_modeled_bacteria_landscape.csv", row.names = F)

#Extract unique phyla and classes
phyla <- gsub("k:Bacteria,p:(\\w+),.*","\\1",hits$V4)
phyla[phyla == "k:Bacteria"] <- "Unknown bacteria"
phyla[phyla == ""] <- "Unknown"
phyla <- gsub("k:Archaea,p:(\\w+),.*","\\1", phyla)
phyla[phyla == "k:Archaea"] <- "Unknown archaea"

class <- gsub("k:Bacteria,p:\\w+,c:(\\w+).*","\\1",hits$V4)
class <- gsub("k:Bacteria,p:\\w+,c:\\[(\\w+).*","\\1", class)
table(class)

#Pie chart
cnvrg_modeled <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG")

head(names(cnvrg_modeled))
tail(names(cnvrg_modeled))
#Subtract 6 from the dimensions for number of taxa present with more than 5 reads.
dim(cnvrg_modeled) - 6
#2360 taxa

#how many were started with:
dat <- read.csv("./processedData/otuTables/smallmem97_16s_for_modeling", 
                fill = T, header = T, stringsAsFactors = F)

dat <- dat[rowSums(dat[,3:length(dat)]) > 0,]
dim(dat) #5532

#how many were modeled
dim(hits) #65

#how many were successfully modeled
results <- read.csv("processedData/all16s.csv")
results <- results[results$taxon != "taxon",]
good <- results[as.numeric(as.character(results$rsq_nested_resampling)) > 0.01,]

# matrix.nrow...1..ncol...1.    taxon rsq_nested_resampling mse_nested_resampling        rsq_reduced          mse_reduced
# 49                       <NA> Zotu2194    0.0513125639367736   0.00108083277982106 0.0158130872822259  0.00110755653284422
# 63                       <NA> Zotu3162    0.0553898902773381   0.00041550342629368  0.029629783810302 0.000415425584516515

########################################
# 2. ITS taxa modeled across all hosts #
########################################

rm(list=ls())

tax <- read.csv("./processedData/ITS_taxa_to_model_via_randomforest.csv")
head(tax)
dim(tax)
sintax <- read.table("./processedData/smallmem97_ITS.sintax", fill = T)

hits <- sintax[sintax$V1 %in% tax$taxa_its,]
write.csv(hits, file = "./processedData/taxonomy_modeled_fungi_landscape.csv", row.names = F)

#Extract unique phyla and classes
phyla <- gsub("d:Fungi,p:(\\w+),.*","\\1",hits$V4)
phyla[phyla == "d:Fungi"] <- "Unknown bacteria"
phyla[phyla == ""] <- "Unknown"

class <- gsub("d:Fungi,p:\\w+,c:(\\w+).*","\\1",hits$V4)
class <- gsub("d:Fungi,p:\\w+,c:\\[(\\w+).*","\\1", class)
table(class)

#Pie chart
cnvrg_modeled <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG")
#how many were started with:
dat <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling", 
                fill = T, header = T, stringsAsFactors = F)

dat <- dat[rowSums(dat[,3:length(dat)]) > 0,]
total <- dim(dat)[1] #3908

#convert data to qualitative data, 1 if present 0 otherwise
dat[,3:length(dat)][dat[,3:length(dat)] > 0] <- 1
dim(dat) #still good
dim(dat[rowSums(dat[,3:length(dat)]) >=100,])[1]


hyperrare <- dim(dat[rowSums(dat[,3:length(dat)]) < 5,])[1]
uncommon <- dim(dat[rowSums(dat[,3:length(dat)]) >= 5 & rowSums(dat[,3:length(dat)]) <= 100 ,])[1]

#how many were successfully modeled
results <- read.csv("processedData/allITS.csv")
results <- results[results$taxon != "taxon",]
#how many were modeled
prevalent <- dim(results)[1] #155

good <- results[as.numeric(as.character(results$rsq_nested_resampling)) >= 0.01,]
predictable <- dim(good)[1]
# 
# matrix.nrow...1..ncol...1.     taxon rsq_nested_resampling mse_nested_resampling          rsq_reduced          mse_reduced
# 21                        <NA>  Zotu1139    0.0413599266381059  6.75687386339465e-06   0.0336250115035115 6.78626707915818e-06
# 35                        <NA>  Zotu1215    0.0368143961339491  1.52955105633748e-05   0.0389912633320888 1.51821825588824e-05
# 55                        <NA> Zotu13140    0.0400502395781242  2.27185373199585e-05  0.00737833068468224 2.24057787841616e-05
# 63                        <NA>  Zotu1340    0.0202105042548996  6.96683174325049e-06    0.014205783192394 6.98338769136251e-06
# 99                        <NA>  Zotu2170    0.0293904847289709   5.8112568522004e-05   0.0404109951457094   5.821917421194e-05
# 131                       <NA>   Zotu329    0.0190225519617162  0.000112312215484468   0.0140571306579045 0.000113668645360282
# 173                       <NA>   Zotu514    0.0539279368002847  5.38460355544884e-05 -0.00423572157557428 5.45495775460423e-05
# 191                       <NA>   Zotu551    0.0745011372353048  3.87301603929449e-05   0.0206415520229568 3.88967782232081e-05
# 197                       <NA>  Zotu5623    0.0866547522474002  2.45225169567343e-05   0.0567346081167543 2.48877820408962e-05
# 213                       <NA>  Zotu5782    0.0356819006587422   0.00011999925067452  -0.0089580603304666 0.000118938306478113
# 245                       <NA>   Zotu677    0.0117982705628421  6.24835105636682e-05  -0.0317356329509936 6.42618405403607e-05
# 249                       <NA>   Zotu687     0.014547665225141  3.06664954622895e-05   0.0121585756131102 2.97144507575842e-05
# 251                       <NA>    Zotu68    0.0123583072345209   0.00145915843264373  0.00717069539876637  0.00146972916678503
# 265                       <NA>   Zotu723     0.010956586259941  2.81025870674192e-05  0.00409466365142658 2.85988889398356e-05
# 271                       <NA>   Zotu735    0.0915126086282025  2.66564625283471e-05   0.0675183178582238 2.63883003524444e-05
# 283                       <NA>    Zotu80    0.0423577505185719  0.000908001721868931   0.0333696809722826 0.000938726184502033
# 285                       <NA>  Zotu8197    0.0301665066081335  2.59388157227417e-05   0.0267335466360699 2.56429430600524e-05

#subtracting stuff so that the total of the pie chart is accurate, things not double counted
slice1 <- hyperrare
slice2 <- uncommon
slice3 <- prevalent - predictable
slice4 <- predictable

slices <- c(slice1, slice2, slice3, slice4)

#QC
sum(slices) == total

lbls <- c("Very rare OTUs  (< 5 counts)",
          "OTUs with 5—29 counts",
          "OTUs with 30+ counts",
          "Prevalent OTUs (in 100+ samples)",
          "Predictable OTUs (R2 > 1%)")

pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels 

pie(slices, labels = lbls, main="")

########################################
# 1. 16s taxa modeled across all hosts #
########################################
