#Determine the extent to which coligos occur where they are not supposed to. 

dat <- read.table("~/Desktop/teton/coligoTable", stringsAsFactors = F, header = T)
dat[1:5,1:5]
dim(dat)

names(dat) <- gsub("_[ABCDEFGH](\\d+)$","", names(dat))

key <- read.csv("~/Desktop/teton/NovaSeq2_DemuxJH.csv", stringsAsFactors = F)

#names(dat) <- gsub("X16S_(.*)","\\1",names(dat))

key$forward_barcode <- toupper(key$forward_barcode)
key$reverse_barcode <- toupper(key$reverse_barcode)

key$combo <- paste(key$locus, key$forward_barcode, key$reverse_barcode, key$samplename, sep = "_")

key$combo[grep("^16S", key$combo)] <- paste("X", key$combo[grep("^16S", key$combo)], sep = "")

#converting names to well position, for ease
for(i in 1:length(key$combo)){
  names(dat)[which(key$combo[i] == names(dat))] <- key$wellposition[i]
}

dat$OTUID <- gsub("Coligo_", "",dat$OTUID)

dat <- dat[dat$OTUID != "ISD",]

percent_target <- NA

for(i in 2:length(names(dat))){
  if(dat[dat$OTUID == names(dat)[i],i] > 15){
  percent_target[i] <- dat[dat$OTUID == names(dat)[i],i] / sum(dat[,i])
  }
}
mean(na.omit(percent_target))
summary(percent_target)
table(percent_target > 0.90)

dat[,1:15]

# IT appears that cases of multiple contamination are pretty high, but the proportion of overall reads 
#that are contaminants are very low. When contamination happened it was normally one well being dragged through
#a whole row


#FOR ITS
#rm(list=ls())
datITS <- read.table("coligoTableITS", stringsAsFactors = F, header = T)
datITS[1:5,1:5]
dim(datITS)

key <- read.csv("novaseq3_demux.csv", stringsAsFactors = F)

names(datITS) <- gsub("ITS_(.*)","\\1",names(datITS))

key$forward_barcode <- toupper(key$forward_barcode)
key$reverse_barcode <- toupper(key$reverse_barcode)

key$combo <- paste(key$forward_barcode, key$reverse_barcode, key$samplename, sep = "_")

for(i in 1:length(key$combo)){
  names(datITS)[which(key$combo[i] == names(datITS))] <- key$wellposition[i]
}

datITS$OTUID <- gsub("Coligo_", "",datITS$OTUID)

datITS <- datITS[datITS$OTUID != "ISD",]

percent_target <- NA

for(i in 2:length(names(datITS))){
  if(datITS[datITS$OTUID == names(datITS)[i],i] > 15){
    percent_target[i] <- datITS[datITS$OTUID == names(datITS)[i],i] / sum(datITS[,i])
  }
}
mean(na.omit(percent_target))
summary(percent_target)
table(percent_target > 0.90)

datITS[,1:15]

#way better performance for ITS
dim(datITS)

#why are there more samples for 16s than ITS? There should be the same amount. 

table(key$locus)
dim(datITS)
#2960 for ITS
#4095 for 16s
#4224 for both is expected

#1.Identify things that didn't show up for ITS
#2. Figure out the coligo for some of those things
#3. spot grep for some coligos in the actual data.
#4. figure out what to do next. 
#Answer: turns out there is a T or A in the front of the coligos, due to faulty primer trimming.
#Need to figure out a better way to do a match. Maybe cut to low complexity data, then remove the primer, then take the
#next 13 or 14 bases or something and loop through looking for specific coligos and do an error correction?
#maybe simply do a two step search-exact where we first lop of the primer and do search_exact then lop off another 
#base and do it again. Super clunky, but would probably work. Obviously would miss cases where the primer was too short. 

dat <- read.table("coligoTable16s", stringsAsFactors = F, header = T)
names(dat) <- gsub("X16S_(.*)","\\1",names(dat))

datITS <- read.table("coligoTableITS", stringsAsFactors = F, header = T)
names(datITS) <- gsub("ITS_(.*)","\\1",names(datITS))
names(datITS) <- gsub("ITS","16S",names(datITS))

table(names(datITS) %in% names(dat))
table(names(dat) %in% names(datITS))

only16 <- names(dat)[which(!(names(dat) %in% names(datITS)))]


key16 <- key[key$combo %in% only16,]
table(key16$midplate)

