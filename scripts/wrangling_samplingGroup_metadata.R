#This script condenses the metadata by sampmling group. 
#this makes it easier to assign colors /subset the data by group during 
#other analyses

rm(list=ls())
#sadly the original treatments metadata file was at least partly made by hand
#during field work, so there is no script to call before this one.
dat <- read.csv("processedData/treatments_metadata.csv", stringsAsFactors = F)
meta <- read.csv("processedData/metadata_2.csv", stringsAsFactors = F)
meta <- meta[meta$substrate == "plant",]
#add compartment column
dat$compartment <- "EP"
dat$compartment[grep("EN", dat$x)] <- "EN"
table(dat$compartment)

dat$plant <- gsub("(\\d+_\\d+_\\d+)E.*", "\\1", dat$x)

dat$taxon <- NA
dat$location <- NA
dat$lifehistory <- NA
for(i in unique(dat$plant)){
  if(length(unique(meta$taxon.y[grep(i, meta$region_site_plant)])) > 1){
    print(i)
  }
  dat$taxon[dat$plant == i] <-
    unique(meta$taxon.y[grep(i, meta$region_site_plant)])
  dat$location[dat$plant == i] <-
    unique(meta$region_site[grep(i, meta$region_site_plant)])
  dat$lifehistory[dat$plant == i] <-
    unique(meta$habit[grep(i, meta$region_site_plant)])
}
head(dat)

dat <- dplyr::distinct(dat[,2:length(dat)])

for(i in unique(dat$plant)){
  if(length(unique(meta$taxon.x[grep(i, meta$region_site_plant)])) > 1){
    print(i)
  }
dat$lifehistory[dat$plant == i] <-
  unique(meta$habit[grep(i, meta$region_site_plant)])
}

dat$lifehistory[dat$taxon == "Picea engelmannii"] <- "tree"
dat$lifehistory[dat$taxon == "Juniperus communis"] <- "shrub"
dat$lifehistory[dat$taxon == "Unknown fir"] <- "tree"
dat$lifehistory[dat$taxon == "Juncus sp."] <- "graminoid"
dat$lifehistory[dat$taxon == "Osmorhiza depauperata"] <- "forb"
dat$lifehistory[dat$taxon == "Pinus contorta"] <- "tree"
dat$lifehistory[dat$taxon == "Paxistima myrsinites"] <- "shrub"
dat$lifehistory[dat$taxon == "Symphoricarpos albus"] <- "shrub"
dat$lifehistory[dat$taxon == "Vaccinium membranaceum 2_1_6"] <- "shrub"
dat$taxon[dat$taxon == "Vaccinium membranaceum 2_1_6"] <- "Vaccinium membranaceum"
dat$lifehistory[dat$taxon == "Potentilla diversifolia var. diversifolia"] <- "forb"
dat$lifehistory[dat$taxon == "Trifolium parryi var. montanense"] <- "forb"
dat$lifehistory[dat$taxon == "Astragalus kentrophyta cf. "] <- "forb"
dat$lifehistory[dat$taxon == "Sedum lanceolatum"] <- "forb"
dat$lifehistory[dat$taxon == "Juncus balticus "] <- "graminoid"
dat$taxon[dat$taxon == "Artemisia tridentate"] <- "Artemisia tridentata"
dat$lifehistory[dat$taxon == "Artemisia tridentate"] <- "shrub"
dat$lifehistory[dat$taxon == "Frasera speciosa"] <- "forb"
dat$lifehistory[dat$taxon == "Antennaria media"] <- "forb"

names(dat)[1] <- "compartment"

#Reorder to match how we modeled it. Load the order that we need
#Check both loci were modeled in same order

treats16s <- read.csv("./processedData/treatments_for_modeling_16s.csv", 
                      stringsAsFactors = F)
treats16s$x <- gsub("(\\d+_\\d+_\\d+E[NP])_.*", "\\1", treats16s$x)

treatsITS <- read.csv("./processedData/treatments_for_modeling_ITS.csv", 
                      stringsAsFactors = F)
treatsITS$x <- gsub("(\\d+_\\d+_\\d+E[NP])_.*", "\\1", treatsITS$x)

table(treats16s == treatsITS) #yay

dat$treatmentClass <- paste(dat$plant, dat$compartment, sep = "")

dat <- dat[match(treats16s$x, dat$treatmentClass),]

table(treats16s$x == dat$treatmentClass)

#####################
#Bring in the neighboring flora diversity
flora <- read.csv(file = "./processedData/shannons_siteLevel_veg.csv",
                  stringsAsFactors = F, header = T)

dat <- read.csv("processedData/treatments_metadata.csv", stringsAsFactors = F)

newdat <- merge(dat, flora, by.y = "newdat.dat.siteLabel",
      by.x = "location")


names(newdat)[length(newdat)] <- "shannons_flora"
write.csv(newdat, file = "./processedData/treatments_metadata.csv", row.names = F)

#Densitometer
denso <- read.csv(file = "./processedData/densitometer_mean_vs_site.csv",
                  stringsAsFactors = F, header = T)

dat <- read.csv("processedData/treatments_metadata.csv", stringsAsFactors = F)

newdat <- merge(denso, dat, 
                by.y = "location",
                by.x = "dat.siteLabel")
write.csv(newdat, file = "./processedData/treatments_metadata.csv", row.names = F)





