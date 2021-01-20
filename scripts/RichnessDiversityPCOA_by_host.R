#Calculate richness and diversity by host. Do an ordination by host.
#Samples are lumped within a site because the comparison of interest is site and host
#and OTU clustering threshold. 

#We will compare shifts in a site across OTU thresholds, so sample number
#will not confound the analysis.
rm(list=ls())
dat <- read.table("processedData/otuTables//its90_otuTable", header = T, stringsAsFactors = F)
#First, check that the dimensions of the table show the same number of samples 
#for each dataset.
if(dim(dat)[1] != 3857){
  print("Warning: dimensions of data do not appear correct")
}

mdat <- read.csv("./processedData/metadata.csv", header = T, stringsAsFactors = F)
mdat <- mdat[mdat$locus == "ITS",]
mdat$combo <- paste(mdat$locus, mdat$forward_barcode, mdat$reverse_barcode, mdat$plant, sep = "_")
dim(mdat)
dim(dat)

#edit the header of the otu table
ENvsEP <- gsub("([ITS16]*_\\w+_\\w+_\\d+_\\d+_\\d+_\\d+)_([ENP]*).*", "\\2", names(dat))
ENvsEP <- gsub("([ITS16]*_\\w+_\\w+_\\d+_\\d+_\\d+_\\d+)([ENPepn]*).*", "\\2",ENvsEP)
ENvsEP <- toupper(ENvsEP)

names(dat) <- gsub("([ITS16]*_\\w+_\\w+_\\d+_\\d+_\\d+_\\d+).*", "\\1", names(dat))
dat <- dat[, names(dat) %in% mdat$combo ]
mdat <- mdat[mdat$combo %in% names(dat),]
ENvsEP <- ENvsEP[order(mdat$combo)] #doublecheck all this reorganizing

mdat <- mdat[order(mdat$combo),]
dat <- dat[,match(mdat$combo, names(dat))]

if(any(names(table(names(dat) == mdat$combo))=="FALSE")){
  print("Warning: metadata and seq data aren't matching")
}

#Estimate richness using Breakaway
#Estimating richness for microbial assemblages is very difficult and frankly 
#I think it is perhaps nonsensical when limited to single-locus sequencing data, but everyone wants to do it....so

library(breakaway)
#To confirm input format check out:
#data(toy_otu_table)
#data(toy_metadata)

frequencytablelist <- build_frequency_count_tables(dat)

#the clunky error handling is needed due to issues with uniroot, a function that gets called by breakaway
rich <- list()
for(i in 1:length(frequencytablelist)){
  if(sum(frequencytablelist[[i]][,1]) > 10 & length(frequencytablelist[[i]][,1]) > 3){ #omit samples without much info
    possibleError <- tryCatch(
      breakaway(frequencytablelist[[i]], cutoff = 100000),
      error=function(e) e)
    if(inherits(possibleError, "error")){
      rich[[i]] <- NA
      next
    }
    rich[[i]] <- breakaway(frequencytablelist[[i]], cutoff = 100000)
  }
}

length(rich)

richest <- NA
for(i in 1:length(rich)){
 if(is.null(rich[[i]])){
   richest[i] <- NA
  }else if( is.na(rich[[i]])){
    richest[i] <- NA
  }else{
    richest[i] <- rich[[i]]$estimate
  }
}
boxplot(log(richest[ENvsEP %in% c("EN","EP")])~ENvsEP[ENvsEP %in% c("EN","EP")])
aggregate(richest[ENvsEP %in% c("EN","EP")]~ENvsEP[ENvsEP %in% c("EN","EP")], FUN = mean)

#Calculate diversity
#Output diversity as a table


#Output ordination of data

