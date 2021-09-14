#two part script one part for each locus
rm(list = ls())
dat <- read.csv("./processedData/multi_ests_16s.csv",
                stringsAsFactors = F)
dat[1:3,1:3]
dim(dat)

metadat <- read.csv("./processedData/metadata_2.csv",
                    stringsAsFactors = F)
table(is.na(metadat$habit))

dim(metadat)

metadat[metadat$plant == "1_1_1_1",] #still good

#CHANGE LOCUS AS NEEDED
metadat <- metadat[metadat$locus == "16S",]

dat$sample <- gsub("X", "",dat$sample)
dim(metadat)

#Need to get metadata into the same order as the modeled data
metadat$sample <- gsub("_16S_[1,2]_\\w+$", "",metadat$sample)
metadat$sample <- gsub("_16S$", "",metadat$sample)
metadat$sample <- gsub("dupe_", "",metadat$sample)
metadat$sample <- gsub("_rewash$", "",metadat$sample)
metadat$sample <- gsub("(E[NP]).*$", "\\1",metadat$sample)
metadat$sample <- gsub("(\\d+)(E[NP])$", "\\1_\\2",metadat$sample)

metadat$samplename[grep("7_1_1_4", metadat$samplename)]

#Remove soil
metadat <- metadat[metadat$substrate == "plant",]

metadat[metadat$plant == "1_1_1_1",20:35] #still good

#There is no need to keep all PCR replicates from the metadata, since they have the 
#same metadata. Therefore, we will go through and keep just one row for each 
#unique sample identifier

metadat_reduced <- data.frame(matrix(nrow = length(unique(metadat$sample)),
                                     ncol = length(metadat)))
k <- 1
for(i in unique(metadat$sample)){
  metadat_reduced[k,] <- metadat[which(metadat$sample == i)[1],]
  k <- k + 1
}
names(metadat_reduced) <- names(metadat)

metadat_reduced[metadat_reduced$plant == "1_1_1_1",20:35] #still good

setdiff(metadat_reduced$sample, dat$sample)

#I spot checked some of these missing things and I removed them during wrangling becaaues
#They were cross-contaminated.
setdiff(dat$sample, metadat_reduced$sample) #nada
dim(dat)
dim(metadat_reduced)

samplemeta <- read.csv("./processedData/treatments_metadata.csv",
                       stringsAsFactors = F)

#Prior to merge, we are good
metadat_reduced[metadat_reduced$plant == "1_1_1_1",] #
samplemeta[samplemeta$plant == "1_1_1",]

names(samplemeta)[names(samplemeta) %in% names(metadat_reduced)]
samplemeta <- samplemeta[,names(samplemeta) !="compartment"]
samplemeta <- samplemeta[!duplicated(samplemeta$plant),]

metadat_reduced2 <- merge(metadat_reduced, 
                         samplemeta, 
                         by.x = "region_site_plant",
                         by.y = "plant")
dim(metadat_reduced2)
dim(metadat_reduced)

metadat_reduced2[metadat_reduced2$plant == "1_1_1_1",1:35] #still good

metadat_reduced2$julianDate <- format(as.Date(gsub("_","-",metadat_reduced2$date),
                                             format = "%m-%d-%Y"), "%j")

#Cancel the suspect thickness data
metadat_reduced2$thickness[metadat_reduced2$thicknessReliable == "N"] <- NA

#bring in abiotic data
abiotic_data <- read.csv("./raw_data_notIncluding_sequenceData/siteData.csv", stringsAsFactors = F)

metadat_reduced3 <-merge(metadat_reduced2, 
                        abiotic_data[,c(1,20:25)], 
                        by.x = "region_site",
                        by.y = "label", 
                        all.x = T)

dim(metadat_reduced2)
dim(metadat_reduced3)
metadat_reduced3[metadat_reduced3$plant == "1_1_1_1",] #still good

#make metadat and dat have same samples
metadat_reduced3 <- metadat_reduced3[metadat_reduced3$sample %in% dat$sample,]
#reorder so that the metadat matches the modeled dat
metadat_reduced3 <- metadat_reduced3[match(dat$sample,
                                          metadat_reduced3$sample),]

table(metadat_reduced3$sample == dat$sample)

metadat_reduced3$plant_vol <- metadat_reduced3$plantMeasurements__height*
  metadat_reduced3$plantMeasurements__width2*metadat_reduced3$width

metadat_reduced3$habit[metadat_reduced3$taxon_final == "Picea engelmannii"] <-  "tree"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Juniperus communis"] <-  "shrub"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Unknown fir"] <-  "tree"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Juncus sp."] <-  "graminoid"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Juncus balticus"] <-  "graminoid"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Osmorhiza depauperata"] <-  "forb"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Paxistima myrsinites"] <-  "shrub"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Astragalus kentrophyta cf."] <-  "forb"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Antennaria media"] <-  "forb"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Frasera speciosa"] <-  "forb"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Trifolium parryi var. montanense"] <-  "forb"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Potentilla diversifolia"] <-  "forb"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Sedum lanceolatum"] <-  "forb"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Symphoricarpos albus"] <-  "shrub"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Vaccinium membranaceum 2_1_6"] <-  "shrub"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Artemisia tridentata"] <-  "shrub"
metadat_reduced3$habit[metadat_reduced3$taxon_final == "Pinus contorta"] <-  "tree"
metadat_reduced3$taxon_final[metadat_reduced3$taxon_final == "Vaccinium membranaceum 2_1_6"] <-"Vaccinium membranaceum"
metadat_reduced3$taxon_final[metadat_reduced3$taxon_final == "Astragalus kentrophyta cf. "] <-  "Astragalus kentrophyta"

table(metadat_reduced3$habit)
table(metadat_reduced3$taxon_final)
table(is.na(metadat_reduced3$taxon_final))
table(is.na(metadat_reduced3$habit))
table(is.na(metadat_reduced3$compartment))

#Convert unknown hosts to known hosts. 
#This should not matter for modeling, since we modeled these separately and these would have been
#an unusual tree at each site, so shouldn't probably have been lumped with other samples

#EDIT: turns out the unknown stuff matches multiple things that could have been present so the
#cpDNA is not good enough to do a strong ID, at least for some taxa. i didn't look up each one 
#as it is a pain, just did a couple

#metadat_reduced3$sample[metadat_reduced3$taxon_final == "Unknown fir"]
#metadat_reduced3$taxon_final[grep("1_3_1_\\d+_E[NP]",metadat_reduced3$sample)] <- 

#Try and fill in the missing compartments
metadat_reduced3$compartment[is.na(metadat_reduced3$compartment)]  <- gsub("[0-9_]*(E[NP])[0-9_ITS]*","\\1",
      metadat_reduced3$samplename[is.na(metadat_reduced3$compartment)])

table(is.na(metadat_reduced3$compartment))

#add densitomter data
denso <- read.csv("processedData/densitometer_mean_vs_site.csv")
denso$dat.siteLabel <- gsub("([1-9]_[1-9]).*","\\1", denso$dat.siteLabel)
denso$dat.siteLabel <- gsub(".*([1-9]_[1-9])","\\1", denso$dat.siteLabel)
dim(metadat_reduced3)

metadat_reduced3 <-
  merge(metadat_reduced3, denso, by.x = "region_site", by.y = "dat.siteLabel", all.x = T)

#Add floral diversity data
veg <- read.csv("processedData/shannons_siteLevel_veg.csv")
veg$newdat.dat.siteLabel <- gsub("([1-9]_[1-9]).*","\\1", veg$newdat.dat.siteLabel)
veg$newdat.dat.siteLabel <- gsub(".*([1-9]_[1-9])","\\1", veg$newdat.dat.siteLabel)

dim(metadat_reduced3)

metadat_reduced3 <-
  merge(metadat_reduced3, veg, by.x = "region_site", by.y = "newdat.dat.siteLabel", all.x = T)
dim(metadat_reduced3)
#commence to doing analyses!
write.csv(metadat_reduced3, file = "./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv", row.names = F)
write.csv(dat, file = "./processedData/16sp_estimates_wrangled_for_post_modeling_analysis.csv", row.names = F)
dat[,3:length(dat)] <- dat[,3:length(dat)] / dat$ISD
write.csv(dat, file = "./processedData/16sp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv", row.names = F)

############
#Do for ITS#
############

rm(list = ls())
  dat <- read.csv("./processedData/multi_ests_its.csv",
                stringsAsFactors = F)
  dat[1:3,1:3]

  metadat <- read.csv("./processedData/metadata_2.csv",
                      stringsAsFactors = F)
  table(is.na(metadat$habit))
  
  dim(metadat)
  
  metadat[metadat$plant == "1_1_1_1",] #still good
  
  #CHANGE LOCUS AS NEEDED
  metadat <- metadat[metadat$locus == "ITS",]
  
  dat$sample <- gsub("X", "",dat$sample)
  dim(metadat)
  
  #Need to get metadata into the same order as the modeled data
  metadat$sample <- gsub("_ITS_[1,2]_\\w+$", "",metadat$sample)
  metadat$sample <- gsub("_iTS$", "",metadat$sample)
  metadat$sample <- gsub("dupe_", "",metadat$sample)
  metadat$sample <- gsub("_rewash$", "",metadat$sample)
  metadat$sample <- gsub("(E[NP]).*$", "\\1",metadat$sample)
  metadat$sample <- gsub("(\\d+)(E[NP])$", "\\1_\\2",metadat$sample)
  
  metadat$samplename[grep("7_1_1_4", metadat$samplename)]
  
  #Remove soil
  metadat <- metadat[metadat$substrate == "plant",]
  
  metadat[metadat$plant == "1_1_1_1",20:35] #still good
  
  #There is no need to keep all PCR replicates from the metadata, since they have the 
  #same metadata. Therefore, we will go through and keep just one row for each 
  #unique sample identifier
  
  metadat_reduced <- data.frame(matrix(nrow = length(unique(metadat$sample)),
                                       ncol = length(metadat)))
  k <- 1
  for(i in unique(metadat$sample)){
    metadat_reduced[k,] <- metadat[which(metadat$sample == i)[1],]
    k <- k + 1
  }
  names(metadat_reduced) <- names(metadat)
  
  metadat_reduced[metadat_reduced$plant == "1_1_1_1",20:35] #still good
  
  setdiff(metadat_reduced$sample, dat$sample)
  
  #I spot checked some of these missing things and I removed them during wrangling becaaues
  #They were cross-contaminated.
  setdiff(dat$sample, metadat_reduced$sample) #nada
  
  
  samplemeta <- read.csv("./processedData/treatments_metadata.csv",
                         stringsAsFactors = F)
  
  #Prior to merge, we are good
  metadat_reduced[metadat_reduced$plant == "1_1_1_1",] #
  samplemeta[samplemeta$plant == "1_1_1",]
  
  names(samplemeta)[names(samplemeta) %in% names(metadat_reduced)]
  samplemeta <- samplemeta[,names(samplemeta) !="compartment"]
  samplemeta <- samplemeta[!duplicated(samplemeta$plant),]
  
  metadat_reduced2 <- merge(metadat_reduced, 
                            samplemeta, 
                            by.x = "region_site_plant",
                            by.y = "plant")
  dim(metadat_reduced2)
  dim(metadat_reduced)
  
  metadat_reduced2[metadat_reduced2$plant == "1_1_1_1",1:35] #still good
  
  metadat_reduced2$julianDate <- format(as.Date(gsub("_","-",metadat_reduced2$date),
                                                format = "%m-%d-%Y"), "%j")
  
  #Cancel the suspect thickness data
  metadat_reduced2$thickness[metadat_reduced2$thicknessReliable == "N"] <- NA
  
  #bring in abiotic data
  abiotic_data <- read.csv("./raw_data_notIncluding_sequenceData/siteData.csv", stringsAsFactors = F)
  
  metadat_reduced3 <-merge(metadat_reduced2, 
                           abiotic_data[,c(1,20:25)], 
                           by.x = "region_site",
                           by.y = "label", 
                           all.x = T)
  
  dim(metadat_reduced2)
  dim(metadat_reduced3)
  metadat_reduced3[metadat_reduced3$plant == "1_1_1_1",] #still good
  
  #make metadat and dat have same samples
  metadat_reduced3 <- metadat_reduced3[metadat_reduced3$sample %in% dat$sample,]
  #reorder so that the metadat matches the modeled dat
  metadat_reduced3 <- metadat_reduced3[match(dat$sample,
                                             metadat_reduced3$sample),]
  
  table(metadat_reduced3$sample == dat$sample)
  
  metadat_reduced3$plant_vol <- metadat_reduced3$plantMeasurements__height*
    metadat_reduced3$plantMeasurements__width2*metadat_reduced3$width
  
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Picea engelmannii"] <-  "tree"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Juniperus communis"] <-  "shrub"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Unknown fir"] <-  "tree"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Juncus sp."] <-  "graminoid"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Juncus balticus"] <-  "graminoid"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Osmorhiza depauperata"] <-  "forb"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Paxistima myrsinites"] <-  "shrub"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Astragalus kentrophyta cf."] <-  "forb"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Antennaria media"] <-  "forb"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Frasera speciosa"] <-  "forb"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Trifolium parryi var. montanense"] <-  "forb"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Potentilla diversifolia"] <-  "forb"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Sedum lanceolatum"] <-  "forb"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Symphoricarpos albus"] <-  "shrub"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Vaccinium membranaceum 2_1_6"] <-  "shrub"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Artemisia tridentata"] <-  "shrub"
  metadat_reduced3$habit[metadat_reduced3$taxon_final == "Pinus contorta"] <-  "tree"
  metadat_reduced3$taxon_final[metadat_reduced3$taxon_final == "Vaccinium membranaceum 2_1_6"] <-"Vaccinium membranaceum"
  metadat_reduced3$taxon_final[metadat_reduced3$taxon_final == "Astragalus kentrophyta cf. "] <-  "Astragalus kentrophyta"
  metadat_reduced3$taxon_final[metadat_reduced3$taxon_final == "Unknown spruce"] <-  "Unknown Spruce"

  table(metadat_reduced3$lifehistory)
  table(metadat_reduced3$taxon_final)
  table(is.na(metadat_reduced3$taxon_final))
  table(is.na(metadat_reduced3$habit))
  table(is.na(metadat_reduced3$compartment))
  
  #Try and fill in the missing compartments
  metadat_reduced3$compartment[is.na(metadat_reduced3$compartment)]  <- gsub("[0-9_]*(E[NP])[0-9_ITS]*","\\1",
                                                                             metadat_reduced3$samplename[is.na(metadat_reduced3$compartment)])
  table(is.na(metadat_reduced3$compartment))
  
  
  #add densitomter data
  denso <- read.csv("processedData/densitometer_mean_vs_site.csv")
  denso$dat.siteLabel <- gsub("([1-9]_[1-9]).*","\\1", denso$dat.siteLabel)
  denso$dat.siteLabel <- gsub(".*([1-9]_[1-9])","\\1", denso$dat.siteLabel)
  dim(metadat_reduced3)
  
  metadat_reduced3 <-
    merge(metadat_reduced3, denso, by.x = "region_site", by.y = "dat.siteLabel", all.x = T)
  
  #Add floral diversity data
  veg <- read.csv("processedData/shannons_siteLevel_veg.csv")
  veg$newdat.dat.siteLabel <- gsub("([1-9]_[1-9]).*","\\1", veg$newdat.dat.siteLabel)
  veg$newdat.dat.siteLabel <- gsub(".*([1-9]_[1-9])","\\1", veg$newdat.dat.siteLabel)
  
  dim(metadat_reduced3)
  
  metadat_reduced3 <-
    merge(metadat_reduced3, veg, by.x = "region_site", by.y = "newdat.dat.siteLabel", all.x = T)
  
#commence to doing analyses!
write.csv(metadat_reduced3, file = "./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv", row.names = F)
write.csv(dat, file = "./processedData/ITSp_estimates_wrangled_for_post_modeling_analysis.csv", row.names = F)
dat[,3:length(dat)] <- dat[,3:length(dat)] / dat$ISD
write.csv(dat, file = "./processedData/ITSp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv", row.names = F)

