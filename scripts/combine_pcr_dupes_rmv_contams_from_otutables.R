#Order of operations in this script is:
#Fix the few plates were I added too much ISD
#Remove cross-contaminated samples
#Examine rewwashed samples
#remove duplicate samples
#remove soil data
#combine pcr dupes
#remove contaminants as determined from blanks

rm(list=ls())

#Note that smallmem* files are clustered ESVs
dat <- read.csv("./processedData/otuTables/16s_97_smallmem_otuTableCLEAN",
                fill = T, header = T, stringsAsFactors = F)
dat <- read.csv("./processedData/otuTables/its97smallmem_otuTableCLEAN",
               fill = T, header = T, stringsAsFactors = F)

#For ease I check for locus here and then wrap code below such that 
#the correct version runs for each locus
if(any(grepl("rna16S", names(dat))==T)){
  locus <- "16s"
}else{
  locus <- "its"
}

#Seeing if the ISD and so on are in the otu table.
#Counts for these were determined within the cleanOTUtable.R program
tail(dat$OTUID)

# Bring in metadata to determine which samples are soil data
metadat <- read.csv("./processedData/metadata.csv", stringsAsFactors = F)

metadat$mid <- paste(metadat$forward_barcode, 
                     metadat$reverse_barcode,
                     sep = "_")
metadat$mid <- paste(metadat$locus,
                     metadat$mid,
                     sep = "")

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

#First remove duds and ISD
dat_rewash <- dat[!(dat$OTUID %in% c("duds", "ISD", "mtDNA", "cpDNA")),]

rewashes <- names(dat_rewash)[grep("rewash", names(dat_rewash))]
length(rewashes) #84 rewashes for ITS; 82 for 16s

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
boxplot(log(rewashTest$ens),
        log(rewashTest$eps),
        log(rewashTest$rewashes),
        names = c("EN","EP","rewash"))
dev.off()

#Unexpected that rewashing led to more reads. 
#Could be that the second wash was more of a good rinse and got more off?
rm(dat_rewash)

###################################
# Examine blanks#
###################################
#these should have way fewer reads in them then non blanks. 


#####################################################################
#get rid of the dupes, make a name vector for removing John's stuff.# 
#####################################################################
dat <- dat[,-grep("dupe", names(dat))]
dat <- dat[,-grep("triplicate", names(dat))]

#Clean up names
names(dat) <- gsub("ep", "_EP", names(dat))
names(dat) <- gsub("en", "_EN", names(dat))
names(dat) <- gsub("^rna","", names(dat))
names(dat) <- gsub("\\d+","", names(dat))
#names(dat) <- gsub("ITS","", names(dat))
names(dat) <- gsub("[Ee]","", names(dat))
names(dat) <- gsub("[Pp]","", names(dat))
names(dat) <- gsub("[Nn]","", names(dat))
names(dat) <- gsub("primary","", names(dat))
names(dat) <- gsub("primaryCHECKpostseq","", names(dat))
names(dat) <- gsub("rimary","", names(dat))

names(dat) <- gsub("_([ATCG]*_[ATCG]*)_.*","\\1", names(dat))

#A few things missing. These things are from Calder's projects, so I will remove them here

if(locus == "16s"){
  #sanity check
  print(which(!(paste("16",names(dat), sep = "") %in% metadat$mid)))
  dat <- dat[,which(paste("16",names(dat), sep = "") %in% metadat$mid)]
  metadat <- metadat[metadat$mid %in% paste("16",names(dat), sep = ""),]
  names(dat) <- paste("16",names(dat), sep = "")
}else if(locus == "its"){
  #For ITS
  print(which(!( names(dat) %in% metadat$mid)))
  dat <- dat[,which(names(dat) %in% metadat$mid)]
  metadat <- metadat[metadat$mid %in% names(dat),]
  }

#Now cut stuff out of metadat so both dat and metadat$mid have the same stuff
#Depending upon the OTU table (e.g., which locus it is for) this will remove
#about half the data.

#sanity check
table(gsub("\\.\\d+","",names(dat)) %in% metadat$mid)
setdiff(gsub("\\.\\d+","",names(dat)) , metadat$mid)
setdiff( metadat$mid, gsub("\\.\\d+","",names(dat)))
dim(dat)
dim(metadat)

#rearrange
metadat <- metadat[metadat$mid %in% gsub("\\.\\d+","",names(dat)),]
# dim(metadat)
# #There are duplicated mids. Why is that?
# metadat[(duplicated(metadat$mid)),]
# #test case
# metadat[grep("ITSTTACGCTCAA_TGAGCTGCCA", metadat$mid),]
# #A: BECAUSE I SEQUENCED on several libraries. 


dat2 <- dat[,match(metadat$mid, 
                   gsub("\\.\\d+","",names(dat)))]

#sanity check to make sure rearrangement worked
table(gsub("\\.\\d+","",names(dat2)) == metadat$mid)

#Remove soil samples (remove 100 samples)
dat2 <- dat2[,metadat$substrate != "soil"]
metadat <- metadat[metadat$substrate != "soil",]
table(gsub("\\.\\d+","",names(dat2)) == metadat$mid)

#combine PCR duplicates
#Nenames is in the same order as samplename in the metadata and the names of dat2
newnames <- gsub("([0-9_]+[ENP]+).*","\\1",metadat$samplename)

newdat <- data.frame(matrix(nrow = dim(dat2)[1]))
k <- 1
l <- 1
checkthese <- NA
for_names <- vector()
for(i in unique(newnames)){
  if(i != "OTUID"){
    #Make sure there are more than one row to add
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
otus[which(perc_contam > 0.01)]
contams[which(perc_contam > 0.01)]

#percentage of total reads that are likely to be contaminants and thus have been removed
sum(contams[which(perc_contam > 0.01)]) /
  sum(noncontams + contams)

#add otus back in to newdat
newdat <- data.frame(otus, newdat)

#Keep duds but remove other contaminants
if(locus == "16s"){
  newdat[-which(perc_contam > 0.01)[1:36],]
}else if (locus == "its"){
  newdat[-which(perc_contam > 0.01),]
}

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
