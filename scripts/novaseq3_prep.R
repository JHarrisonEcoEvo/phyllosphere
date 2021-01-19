
rm(list=ls())

#I originally made the barcode sheet with the R program below, however it used a sheet in an excel spreadsheet
#tht had incorrect midplate to mid sequence pairings. I am remedying that here. 
dat <- read.csv("~/Desktop/novaseq3_demux.csv", header = T, stringsAsFactors = F)
midpairs <- read.csv("~/Downloads/16SBarcodePairs.csv", header = T, stringsAsFactors = F)

for(i in unique(dat$midplate)){
  for(m in unique(dat$wellposition)){
    if(any(i == dat$midplate & m == dat$wellposition)){
      dat$forward_barcode[i == dat$midplate & m == dat$wellposition] <- midpairs$forward_barcode[i == midpairs$PlateCode & m == midpairs$Position]
      dat$reverse_barcode[i == dat$midplate & m == dat$wellposition] <- midpairs$reverse_barcode[i == midpairs$PlateCode & m == midpairs$Position]
    }
  }
}

head(dat)
table(dat$locus)
table(dat$project)
table(dat$wellposition)
table(dat$plate)
table(dat$client_name)
table(dat$substrate)
table(dat$library)
table(dat$midplate)

#dat <- dat[dat$library == "Novaseq3",]

table(duplicated(dat$samplename))
table(duplicated(paste(dat$forward_barcode, dat$reverse_barcode, dat$locus)))

dat[duplicated(paste(dat$forward_barcode, dat$reverse_barcode, dat$locus)),]

#dat <- dat[,-length(dat)]

grep("\\-", dat$samplename)
grep("\\.", dat$samplename)

write.csv(dat, file = "~/Desktop/novaseq3_demux.csv", row.names = F)


#################
# OLD STUFF BLW #
#################


# dat <- read.csv("~/Desktop/improvedmidkey.csv")
# 
# #Doublecheck that the same capitalization and spelling has been used for factors
# table(dat$library)
# table(dat$locus)
# table(dat$EN_vs_EP)
# 
# dat$sample <- gsub("(\\d+)dupe", "\\1_dupe", dat$sample)
# 
# # unique(dat$notes)
# # dat$rewash <- "primary"
# # dat$rewash[grep("rewash",dat$notes)] <- "rewash"
# 
# dat$sample_name_for_demux <- paste(dat$sample, dat$EN_vs_EP, dat$locus, dat$PCR_duplicate_num, dat$rewash, sep = "_")
# 
# #sort by well, then by plate, to get into the same order as Gregg's mid key for orderliness sake
# dat$well_number <- gsub("[A-Za-z](\\d+)","\\1",dat$well)
# dat <- dat[order(as.numeric(dat$well_number)),]
# 
# dat <- dat[order(dat$Mid_plate_name),]
# dat <- dat[order(dat$locus),]
# 
# dat$well <- as.character( dat$well)
# dat$well <- toupper(dat$well)
# write.csv(dat, file = "~/Desktop/improvedmidkey.csv", row.names = F)
# 
# #merge with barcode key
# og_mid_key <- readxl::read_excel("~/Desktop/16S_ITS_BarcodePairs.xlsx")
# 
# #make columns to merge by
# og_mid_key$combo <- paste(og_mid_key$PlateCode, og_mid_key$Position)
# dat$combo <- paste(dat$Mid_plate_name, dat$well)
# 
# dim(dat)
# mergeddat <- merge(dat, og_mid_key, by.x = "combo", by.y = "combo", all.x = T)
# dim(mergeddat)
# 
# #now convert 16S to ITS in og key and merge again
# og_mid_key$PlateCode <- gsub("16S", "ITS", og_mid_key$PlateCode)
# og_mid_key$combo <- paste(og_mid_key$PlateCode, og_mid_key$Position)
# 
# mergeddat <- merge(mergeddat, og_mid_key, by.x = "combo", by.y = "combo", all.x = T)
# dim(mergeddat)
# 
# #write to excel and clean up a bit
# 
# write.csv(mergeddat, file = "mergedMIDkey.csv", row.names = F)
 dat <- read.csv("~/Desktop/mergedMIDkey.csv")

#subset to each library for ease
novaseq2 <- dat[dat$library == "Novaseq2",]
novaseq3 <- dat[dat$library == "Novaseq3",]

###############################################################################
# Check that there aren't any dupes of mids x locus pairings within a library #
###############################################################################

#Paste mids and locus together and make sure there are no dupes.
combinedMids_novaseq2 <- paste(novaseq2$locus, novaseq2$ForwardMid, novaseq2$ReverseMid, sep = "")

#Dupes are NAs
table(duplicated(combinedMids_novaseq2))
combinedMids_novaseq2[duplicated(combinedMids_novaseq2)]

#Do over for Novaseq 3
combinedMids_novaseq3 <- paste(novaseq3$sample, novaseq3$locus, novaseq3$ForwardMid, novaseq3$ReverseMid, sep = "")

table(duplicated(combinedMids_novaseq3))

#########################################
# Check sample names for demux #
#########################################

#Use this to check that we don't have the wrong PCR duplicate assigned to any samples, should all = 192
table(paste(novaseq2$PCR_duplicate_num, novaseq2$plate)) == 192
table(paste(novaseq3$PCR_duplicate_num, novaseq3$plate)) == 192
which(table(novaseq3$Mid_plate_name) > 96)
which(table(novaseq2$Mid_plate_name) > 96)

#find names that are duplicated
table(duplicated(novaseq3$sample_name_for_demux))

#why dupes?
#check original key

ogKey <- readxl::read_excel("~/Downloads/sampleWell_key_extraction.xlsx")
table(duplicated(paste(ogKey$sample, ogKey$en_vs_ep)))
table(duplicated(dat$sample_name_for_demux)) #same numbers pretty much, so were duplicated off the bat

paste(dat$sample, dat$EN_vs_EP)[table(paste(dat$sample, dat$EN_vs_EP)) > 4]

#do they have unique mids? Yes apparently. I guess we knew this already from above
table(table(paste(dat$sample, dat$EN_vs_EP, dat$locus, dat$ForwardMid, dat$ReverseMid)) > 1) #should be 2 bc of two loci

#See if duplicates have ep /en matches, maybe they are mislabeld

dupedNames <- paste(dat$sample_name_for_demux)[duplicated(paste(dat$sample_name_for_demux))]

#For en samples
#This converts EN to EP then searches our names to see if the EP exists. In many cases it didn't
# I will check after sequencing and see if I can figure out what went down with these. For now I am appending CHECK_post_seq to the sample
gsub("EN", "EP",dupedNames[grep("EN",dupedNames)])[
  !(gsub("EN", "EP",dupedNames[grep("EN",dupedNames)]) %in% dat$sample_name_for_demux)]

gsub("EP", "EN",dupedNames[grep("EP",dupedNames)])[
  !(gsub("EP", "EN",dupedNames[grep("EP",dupedNames)]) %in% dat$sample_name_for_demux)]

dat$sample_name_for_demux  <- as.character(dat$sample_name_for_demux)

#dat$sample_name_for_demux[duplicated(dat$sample_name_for_demux )] <-
#  paste(dat$sample_name_for_demux[duplicated(dat$sample_name_for_demux )], "_CHECK_post_seq", sep = "")

#Doublecheck no dupes now
table(duplicated(dat$sample_name_for_demux ))
dat$sample_name_for_demux[duplicated(dat$sample_name_for_demux )]
#just a few, will deal with in Excel

write.csv(dat, "~/Desktop/mergedMIDkey.csv")

#ogKey <- read.csv("~/Downloads/sampleWell_key_extraction - Sheet1.csv")
#table(duplicated(paste(ogKey$sample, ogKey$en_vs_ep, ogKey$notes)))
