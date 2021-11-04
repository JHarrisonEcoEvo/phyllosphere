#Order of operations in this script is:
#Fix the few plates were I added too much ISD
#Remove cross-contaminated samples
#Examine rewwashed samples
#remove duplicate samples
#remove soil data
#combine pcr dupes
#remove contaminants as determined from blanks

rm(list=ls())

# # #Note that smallmem* files are clustered ESVs
#dat <- read.csv("./processedData/otuTables/16s_97_smallmem_otuTableCLEAN",
#                fill = T, header = T, stringsAsFactors = F)
  dat <- read.csv("./processedData/otuTables/its97smallmem_otuTableCLEAN",
                 fill = T, header = T, stringsAsFactors = F)

 
 #For ease I check for locus here and then wrap code below such that 
 #the correct version runs for each locus
 if(any(grepl("rna16S", names(dat))==T)){
   locus <- "16s"
 }else{
   locus <- "its"
 }
 
#remove soil and other extraneous samples
dat <- dat[,-(grep("[OM]$", names(dat)))]
dat <- dat[,names(dat) != "rna16S_TTGGTAAGAA_CTGCAGACCA_3_3_3_10"]
dat <- dat[,names(dat) != "rna16S_TCCTCTTGAA_CTGCAGACCA_3_3_3_10"]
dat <- dat[, -grep("2_1_1_2", names(dat))] # Remove 2_1_1_2 as I don't know is EN and which EP. Got mislabeled. 

if(locus == "16s"){
  dat <- dat[, -grep("p[23]", names(dat))]
}
dat <- dat[, -grep("_19$", names(dat))]

names(dat) <- gsub("(\\d+)E", "\\1_E", names(dat))
#erasing dupes. These should get added together during the combining of duplicates part below
names(dat) <- gsub("dupe", "", names(dat))

#Seeing if the ISD and so on are in the otu table.
#Counts for these were determined within the cleanOTUtable.R program
tail(dat$OTUID)

# Bring in metadata to determine which samples are soil data
metadat <- read.csv("./processedData/metadata.csv", stringsAsFactors = F)

############################
# Fix over-addition of ISD #
############################
#During early lab work I accidentally added too much ISD to some plates. 

isd <- dat[which(dat[,1] == "ISD"),]
# hist(unlist(isd[,2:length(isd)]))

#Extract those samples that were on problematic plates

prob5k <- metadat$combo[metadat$plate %in% c("19", "14", "9")]
prob3.33 <- metadat$combo[metadat$plate %in% c("23")]

#Writing these to disk. That way I can see if outliers in ordinations or something
#are these samples. 

if(locus == "16s"){
  write.csv(c(prob5k, prob3.33), file = "./processedData/16Ssamples_that_had_too_much_ISD_added.csv", row.names = F)
  
  good <- isd[,!(gsub("rna", "", names(dat)) %in% c(prob5k, prob3.33))]
  print(summary(unlist(good[,2:length(good)])))
  
  #compare to bad samples
  print(summary(unlist(isd[,gsub("rna", "", names(dat)) %in% prob5k])))
  print(summary(unlist(isd[,gsub("rna", "", names(dat))  %in% prob3.33])))
  
  #for 16s
  isd[,gsub("rna", "",names(dat)) %in% prob5k] <- isd[,gsub("rna", "",names(dat)) %in% prob5k] / 5000
  isd[,gsub("rna", "",names(dat)) %in% prob3.33] <- isd[,gsub("rna", "",names(dat)) %in% prob3.33] / 3.33
  
}else if(locus == "its"){
  write.csv(c(prob5k, prob3.33), file = "./processedData/ITSsamples_that_had_too_much_ISD_added.csv", row.names = F)
  
  #See what our counts for the ISD were for good samples
  good <- isd[,!(names(dat) %in% c(prob5k, prob3.33))]
  print(summary(unlist(good[,2:length(good)])))
  
  print(summary(unlist(isd[,names(dat) %in% prob5k])))
  print(summary(unlist(isd[,names(dat) %in% prob3.33])))
  #Much higher in bad samples of course
  
  #Divide problematic samples by how much extra ISD I added
  #This should not be compromised by combining PCR replicates
  isd[,names(dat) %in% prob5k] <- isd[,names(dat) %in% prob5k] / 5000
  isd[,names(dat) %in% prob3.33] <- isd[,names(dat) %in% prob3.33] / 3.33
}else{
  "there is a problem with locus"
}

#replace the ISD in the data frame
dat[which(dat[,1] == "ISD"),] <- isd

###################################
#Remove cross contaminated samples #
###################################

#Remove stuff that is contaminated as shown by coligos.
#These files were made using the coligoChecker.R script
contams <- read.csv("processedData/cross_contam_oct9", stringsAsFactors = F)
contams2 <- read.csv("processedData/cross_contam_sep17", stringsAsFactors = F)
contams <- c(contams$x, contams2$x)

if(locus == "16s"){
  contams <- contams[grep("16S", contams)]
  #16s
  contams %in% paste("rna",metadat$combo, sep = "") 
  metadat <- metadat[!(paste("rna",metadat$combo, sep = "") %in% contams),]
}else if(locus == "its"){
  #ITS
  contams <- contams[grep("ITS", contams)]
  
  contams %in% paste("",metadat$combo, sep = "") 
  metadat <- metadat[!(paste("",metadat$combo, sep = "") %in% contams),]
}

#remove contams
dat <- dat[,!(names(dat) %in% contams)]

###################################
# Examine rewashed samples #
###################################

#normalize
for(i in 2:length(dat)){
  dat[,i] <-  dat[,i] / dat[dat$OTUID == "ISD",i]
}

#First remove duds and ISD
dat_rewash <- dat[!(dat$OTUID %in% c("duds", "ISD", "mtDNA", "cpDNA")),]

rewashes <- names(dat_rewash)[grep("rewash", names(dat_rewash))]
length(rewashes) 


#Compare the rewashes to the original samples. 
#First. grep the prefix of the sample name to get all replicates of the sample.
#second, do count sums of the raw data for the OG data and the rewashed data. 
#May need to take the average across the multiple OG replicates for a sample. 
#save the sample name
#Save these sums to vectors, for EN, EP, rewash, and do boxplot. 

rewashTest <- data.frame(matrix(nrow = length(rewashes)))
k <- 1
for(i in rewashes){
  sample <- gsub(".*_(\\d+_\\d+_\\d+_\\d+).*", "\\1", i)
  hits <- dat_rewash[,grep(sample, names(dat_rewash))]
  if(length(grep("[Ee][Nn]", names(hits))) > 1){
    rewashTest$ens[k] <- mean(colSums(hits[,grep("[Ee][Nn]", names(hits))]))
  }else if(any(grepl("[Ee][Nn]", names(hits)))){
    rewashTest$ens[k] <- sum(hits[,grep("[Ee][Nn]", names(hits))])
  }
  if(length(grep("rewash", names(hits))) > 1){
    rewashTest$rewashes[k] <- mean(colSums(hits[,grep("rewash", names(hits))])) 
  }else if(any(grepl("rewash", names(hits)))){ #There will always be a rewash, but using this code for consistency
    rewashTest$rewashes[k] <- sum(hits[,grep("rewash", names(hits))]) 
  }
  #Remove the rewashes so I don't count them as EPs
  hits <- hits[,-grep("rewash", names(hits))]
  
  if(length(grep("[Ee][Pp]", names(hits))) > 1){
    rewashTest$eps[k] <- mean(colSums(hits[,grep("[Ee][Pp]", names(hits))]))
  }else if(any(grepl("[Ee][Pp]", names(hits)))){
    rewashTest$eps[k] <- sum(hits[,grep("[Ee][Pp]", names(hits))])
  }
  rewashTest$sampleName[k] <- sample
  k <- k + 1
}

pdf(width = 6, height = 6, file = paste("./visuals/rewashBoxplot_", locus,".pdf", sep = ""))
boxplot(log10(rewashTest$ens),
        log10(rewashTest$eps),
        log10(rewashTest$rewashes),
        names = c("EN","EP","rewash"),
        outline = F,
        ylab = "log 10 ISD normalized counts")
dev.off()

#Unexpected that rewashing led to more reads. 
#Could be that the second wash was more of a good rinse and got more off?

#Lets do an ordination. 
#we might hypothesize the rewashes to be  in between EP and EN. 
#This is true for 16s but not for ITS. Maybe more fungi grow both in and out of the leaf than do bacteria.

#grep multiple patterns

matches <- unique(grep(paste(rewashTest$sampleName,collapse="|"), 
                        names(dat_rewash), value=TRUE))

forord <- dat_rewash[,names(dat_rewash) %in% matches]
library(vegan)

#Transposing matrix so vegan functions work. 
#It needs taxa as columns
datr <-  decostand(t(forord), method = "hellinger")

dat_e <- dist(datr,
              method = "euclidean")

ord <- cmdscale(dat_e, 
                eig = TRUE, 
                k = 2)

cols <- vector(mode = "character", 
              length = length(names(forord)))
cols[grep("[Ee][Nn]",names(forord))] <- "orange"
cols[grep("[Ee][Pp]",names(forord))] <- "turquoise"
#Then colorize rewashes, thus overwriting the Ep color
cols[grep("rewash",names(forord))] <- "purple"

#Permanova
adonisOut <- adonis2(dat_e ~ cols, 
                     permutations = 999, method="euclidean")

pdf(width = 6.5, height = 6.5, file = paste("./visuals/ordination_rewashes_", locus,".pdf", sep = ""))

plot(ord$points[,1],ord$points[,2],
     pch = 19,
     col = cols,
     ylab = "PCoA 2",
     xlab = "PCoA 1")
legend("topleft",legend = c("EN", "EP", "rewash"),
       pch = 15,
       col = c("orange","turquoise","purple"))

dev.off()

#names(forord)[cols == ""]

###################################
# Examine blanks#
###################################
#these should have way fewer reads in them then non blanks. 

#First remove duds and ISD

blanks <- names(dat_rewash)[grep("blank", names(dat_rewash))]
length(blanks) #22 blanks

#Almost nothing. Excellent.
colSums(dat_rewash[,names(dat_rewash) %in% blanks])
#Lets try with everything but ISD.
colSums(dat[dat$OTUID != "ISD",names(dat) %in% blanks])

rm(dat_rewash)

########################################################
# Get the otu table in the same format as the metadata #
########################################################

#Will use this later
otus <- dat$OTUID

if(locus == "16s"){
  #sanity check
  names(dat)[which(!(gsub("rna", "",names(dat)) %in% metadat$combo))]
  
  dat <- dat[,which(gsub("rna", "",names(dat)) %in% metadat$combo)]
  metadat <- metadat[metadat$combo %in% gsub("rna", "",names(dat)),]
  #rearrange
  dat2 <- dat[,match(metadat$combo, 
                     gsub("rna", "",names(dat)))]
  #sanity check to make sure rearrangement worked
  table(gsub("rna", "",names(dat2)) == metadat$combo)
}else if(locus == "its"){
  #For ITS
  names(dat)[which(!( names(dat) %in% metadat$combo))]
  dat <- dat[,which(names(dat) %in% metadat$combo)]
  metadat <- metadat[metadat$combo %in% names(dat),]
  #rearrange
  dat2 <- dat[,match(metadat$combo,names(dat))]
  #sanity check to make sure rearrangement worked
  table(names(dat2) == metadat$combo)
}

#Remove soil samples (remove 100 samples)
dat2 <- dat2[,metadat$substrate != "soil"]
metadat <- metadat[metadat$substrate != "soil",]

########################
#combine PCR duplicates#
########################

#Nenames is in the same order as samplename in the metadata and the names of dat2
newnames <- gsub("([0-9_]+[ENP]+).*","\\1",metadat$samplename)

newdat <- data.frame(matrix(nrow = dim(dat2)[1]))
k <- 1
l <- 1
checkthese <- NA
for_names <- vector()
for(i in unique(newnames)){
  if(i != "OTUID"){
    #Make sure there is more than one row to add
    if(length(grep(paste("^",i,"$", sep = ""), newnames)) > 1){
      newdat[,k] <- rowSums(dat2[,grep(paste("^",i,"$", sep = ""), newnames)])
      for_names <- c(for_names, i)
      k <- k + 1
    }
  }
}

dim(newdat)
names(newdat) <- for_names

#Extract blanks
contams <- rowSums(newdat[,grep("lank", names(newdat))])
noncontams <- rowSums(newdat[,-grep("lank", names(newdat))])

perc_contam <- contams / (contams + noncontams)

#percentage of total reads that are likely to be contaminants and thus have been removed
sum(contams[which(perc_contam > 0.01)]) /
  sum(noncontams + contams)

#add otus back in to newdat
newdat <- data.frame(otus, newdat)

#Keep duds but remove other contaminants
duds <- newdat[which(as.character(newdat$otus) == "duds"),]
ISD <- newdat[which(as.character(newdat$otus) == "ISD"),]
mtDNA <- newdat[which(as.character(newdat$otus) == "mtDNA"),]
cpDNA <- newdat[which(as.character(newdat$otus) == "plant"),]

newdat <- newdat[-which(perc_contam > 0.01),]
#newdat[1 + length(newdat$otus),] <- duds
if(locus = "16s"){
  newdat[1 + length(newdat$otus),] <- ISD
}
# newdat[1 + length(newdat$otus),] <- mtDNA
# newdat[1 + length(newdat$otus),] <- cpDNA

otus <- newdat[,1]
newdat <- newdat[,-1]
newdat <- newdat[,-grep("lank", names(newdat))]
dim(newdat)

#######################
#organize for modeling
#######################

sample <- gsub("X(\\d+_\\d+_\\d+_\\d+)_.*","\\1", names(newdat))
compartment <- gsub("X\\d+_\\d+_\\d+_\\d+_(ENP)*","\\1", names(newdat))

newdat2 <- newdat[,order(compartment,sample)]
newdat2 <- data.frame(otus, newdat2)

if(locus == "16s"){
  write.csv(newdat2, file = paste("smallmem97_16s", "for_modeling", sep = "_"))
}else if (locus == "its"){
  write.csv(newdat2, file = paste("smallmem97_ITS", "for_modeling", sep = "_"))
}
