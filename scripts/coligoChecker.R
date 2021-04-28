#Determine the extent to which coligos occur where they are not supposed to. 
#Go by library
rm(list=ls())
dat <- read.table("./processedData/coligoTable_oct9", stringsAsFactors = F, header = T)
dat[1:5,1:5]
dim(dat)

#Novaseq 2 is for Sept 17 library and 3 is the Oct library
nms <- names(dat)
#key <- read.csv("./raw_data_notIncluding_sequenceData/demultiplexingKeys/NovaSeq2_DemuxJH.csv",
#                stringsAsFactors = F)
#dim(key)
#key <- key[,names(key) != "key"]
key <- read.csv("./raw_data_notIncluding_sequenceData/demultiplexingKeys/novaseq3_demux.csv",
                   stringsAsFactors = F)

#names(dat) <- gsub("X16S_(.*)","\\1",names(dat))

key$forward_barcode <- toupper(key$forward_barcode)
key$reverse_barcode <- toupper(key$reverse_barcode)

key$combo <- paste(key$locus, key$forward_barcode, key$reverse_barcode,
                    sep = "_")
key$combo[grep("16S", key$combo)] <- paste("rna", key$combo[grep("16S", key$combo)], sep = "")

#Relic from a failed idea
# dupes <- key$combo[duplicated(key$combo)]
# length(dupes)
# dupes <- key[key$combo %in% dupes,]
# 
# out <- NA
# for(i in dupes$combo){
#   x <- dupes$wellposition[dupes$combo == i]
#   out[i] <- x[1] == x[2]
#   if(out[i] == F){
#     print(x)
#   }
# }
# table(out)

#convert names to well position
for(i in 1:length(key$combo)){
  if(length(grep(key$combo[i], names(dat)) == 1)){
    names(dat)[grep(key$combo[i], names(dat))] <- key$wellposition[i]
  }else{
    next
  }
}

dat$OTUID <- gsub("Coligo_", "",dat$OTUID)

dat <- dat[dat$OTUID != "ISD",]

percent_target <- NA

for(i in 2:length(names(dat))){
  if(names(dat)[i] %in% dat$OTUID){
    if(dat[dat$OTUID == names(dat)[i],i] >= 15){
      percent_target[i] <- dat[dat$OTUID == names(dat)[i],i] / sum(dat[,i])
    }
  }
}
mean(na.omit(percent_target))
summary(percent_target)
table(percent_target > 0.80)

#sanity check
dat[,nms == "rna16S_TCTCTCCAA_CTTCAATAA_5_3_2_6_EP_16S_1_primary"]

write.csv(nms[which(percent_target <= 0.5)], 
          file = "./processedData/cross_contam_oct9",
          row.names = F)

################
#Next library ##
################
rm(list=ls())
dat <- read.table("./processedData/coligoTable_sep17", stringsAsFactors = F, header = T)
dat[1:5,1:5]
dim(dat)

#Novaseq 2 is for Sept 17 library and 3 is the Oct library
nms <- names(dat)
key <- read.csv("./raw_data_notIncluding_sequenceData/demultiplexingKeys/NovaSeq2_DemuxJH.csv",
                stringsAsFactors = F)

key$forward_barcode <- toupper(key$forward_barcode)
key$reverse_barcode <- toupper(key$reverse_barcode)

key$combo <- paste(key$locus, key$forward_barcode, key$reverse_barcode,
                   sep = "_")
key$combo[grep("16S", key$combo)] <- paste("rna", key$combo[grep("16S", key$combo)], sep = "")

#Relic from a failed idea
# dupes <- key$combo[duplicated(key$combo)]
# length(dupes)
# dupes <- key[key$combo %in% dupes,]
# 
# out <- NA
# for(i in dupes$combo){
#   x <- dupes$wellposition[dupes$combo == i]
#   out[i] <- x[1] == x[2]
#   if(out[i] == F){
#     print(x)
#   }
# }
# table(out)

#convert names to well position
for(i in 1:length(key$combo)){
  if(length(grep(key$combo[i], names(dat)) == 1)){
    names(dat)[grep(key$combo[i], names(dat))] <- key$wellposition[i]
  }else{
    next
  }
}

dat$OTUID <- gsub("Coligo_", "",dat$OTUID)

dat <- dat[dat$OTUID != "ISD",]

percent_target <- NA

for(i in 2:length(names(dat))){
  if(names(dat)[i] %in% dat$OTUID){
    if(dat[dat$OTUID == names(dat)[i],i] >= 15){
      percent_target[i] <- dat[dat$OTUID == names(dat)[i],i] / sum(dat[,i])
    }
  }
}
mean(na.omit(percent_target))
summary(percent_target)
table(percent_target > 0.80)

#sanity check
dat[,nms == "rna16S_TGAGAAGAA_TTGCGAAC_4_1_6_7ep"]

write.csv(nms[which(percent_target <= 0.5)], 
          file = "./processedData/cross_contam_sep17",
          row.names = F)
